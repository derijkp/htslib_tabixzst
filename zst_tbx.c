/*  tbx.c -- tabix API functions.

    Copyright (C) 2009, 2010, 2012-2015, 2017-2020, 2022-2023 Genome Research Ltd.
    Copyright (C) 2010-2012 Broad Institute.

    Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#define HTS_BUILDING_LIBRARY // Enables HTSLIB_EXPORT, see htslib/hts_defs.h
#include <config.h>

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <errno.h>
#include "htslib/tbx.h"
#include "htslib/bgzf.h"
#include "htslib/hts_endian.h"
#include "hts_internal.h"
#include "zstdtools.h"

#include "htslib/khash.h"
KHASH_DECLARE(s2i, kh_cstr_t, int64_t)

typedef struct {
    int64_t beg, end;
    char *ss, *se;
    int tid;
} tbx_intv_t;

/*
 * Return a virtual file pointer to the current location in the file.
 * No interpetation of the value should be made, other than a subsequent
 * call to tbx_zst_seek can be used to position the file at the same point.
 * Return value is non-negative on success.
 * Returns -1 on error.
 */
#define zst_tbx_tell(fp) ((fp->framefilepos << 16) | ((fp->outbufferpos) & 0xFFFF))

// just copy the unexported functions used (that are still exactly the same) here for now

static inline int zst_get_tid(tbx_t *tbx, const char *ss, int is_add)
{
    khint_t k;
    khash_t(s2i) *d;
    if (tbx->dict == 0) tbx->dict = kh_init(s2i);
    if (!tbx->dict) return -1; // Out of memory
    d = (khash_t(s2i)*)tbx->dict;
    if (is_add) {
        int absent;
        k = kh_put(s2i, d, ss, &absent);
        if (absent < 0) {
            return -1; // Out of memory
        } else if (absent) {
            char *ss_dup = strdup(ss);
            if (ss_dup) {
                kh_key(d, k) = ss_dup;
                kh_val(d, k) = kh_size(d) - 1;
            } else {
                kh_del(s2i, d, k);
                return -1; // Out of memory
            }
        }
    } else k = kh_get(s2i, d, ss);
    return k == kh_end(d)? -1 : kh_val(d, k);
}

int tbx_parse1(const tbx_conf_t *conf, size_t len, char *line, tbx_intv_t *intv);

static inline int zst_get_intv(tbx_t *tbx, kstring_t *str, tbx_intv_t *intv, int is_add)
{
    if (tbx_parse1(&tbx->conf, str->l, str->s, intv) == 0) {
        int c = *intv->se;
        *intv->se = '\0'; intv->tid = zst_get_tid(tbx, intv->ss, is_add); *intv->se = c;
        if (intv->tid < 0) return -2;  // zst_get_tid out of memory
        return (intv->beg >= 0 && intv->end >= 0)? 0 : -1;
    } else {
        char *type = NULL;
        switch (tbx->conf.preset&0xffff)
        {
            case TBX_SAM: type = "TBX_SAM"; break;
            case TBX_VCF: type = "TBX_VCF"; break;
            case TBX_UCSC: type = "TBX_UCSC"; break;
            default: type = "TBX_GENERIC"; break;
        }
        hts_log_error("Failed to parse %s, was wrong -p [type] used?\nThe offending line was: \"%s\"",
            type, str->s);
        return -1;
    }
}

static int zst_tbx_set_meta(tbx_t *tbx)
{
    int i, l = 0, l_nm;
    uint32_t x[7];
    char **name;
    uint8_t *meta;
    khint_t k;
    khash_t(s2i) *d = (khash_t(s2i)*)tbx->dict;

    memcpy(x, &tbx->conf, 24);
    name = (char**)malloc(sizeof(char*) * kh_size(d));
    if (!name) return -1;
    for (k = kh_begin(d), l = 0; k != kh_end(d); ++k) {
        if (!kh_exist(d, k)) continue;
        name[kh_val(d, k)] = (char*)kh_key(d, k);
        l += strlen(kh_key(d, k)) + 1; // +1 to include '\0'
    }
    l_nm = x[6] = l;
    meta = (uint8_t*)malloc(l_nm + 28);
    if (!meta) { free(name); return -1; }
    if (ed_is_big())
        for (i = 0; i < 7; ++i)
            x[i] = ed_swap_4(x[i]);
    memcpy(meta, x, 28);
    for (l = 28, i = 0; i < (int)kh_size(d); ++i) {
        int x = strlen(name[i]) + 1;
        memcpy(meta + l, name[i], x);
        l += x;
    }
    free(name);
    hts_idx_set_meta(tbx->idx, l, meta, 0);
    return 0;
}

// Minimal effort parser to extract reference length out of VCF header line
// This is used only used to adjust the number of levels if necessary,
// so not a major problem if it doesn't always work.
static void zst_adjust_max_ref_len_vcf(const char *str, int64_t *max_ref_len)
{
    const char *ptr;
    int64_t len;
    if (strncmp(str, "##contig", 8) != 0) return;
    ptr = strstr(str + 8, "length");
    if (!ptr) return;
    for (ptr += 6; *ptr == ' ' || *ptr == '='; ptr++) {}
    len = strtoll(ptr, NULL, 10);
    if (*max_ref_len < len) *max_ref_len = len;
}

// Same for sam files
static void zst_adjust_max_ref_len_sam(const char *str, int64_t *max_ref_len)
{
    const char *ptr;
    int64_t len;
    if (strncmp(str, "@SQ", 3) != 0) return;
    ptr = strstr(str + 3, "\tLN:");
    if (!ptr) return;
    ptr += 4;
    len = strtoll(ptr, NULL, 10);
    if (*max_ref_len < len) *max_ref_len = len;
}

// Adjusts number of levels if not big enough.  This can happen for
// files with very large contigs.
static int zst_adjust_n_lvls(int min_shift, int n_lvls, int64_t max_len)
{
    int64_t s = 1LL << (min_shift + n_lvls * 3);
    max_len += 256;
    for (; max_len > s; ++n_lvls, s <<= 3) {}
    return n_lvls;
}



int zstd_tbx_getline(ZSTDres *res, int delim, kstring_t *str) {
	char *writebuffer;
	uint64_t read;
	unsigned int writesize;
	str->l = 0;
	if (res->contentsize == 0) {
		while (res->contentsize == 0) {
			if (feof(res->finput)) break;
			zstd_skipframe(res);
			zstd_readheader(res);
		}
		zstd_readframe(res);
		res->outbufferpos = 0;
	} else if (!res->frameread) {
		zstd_readframe(res);
		res->outbufferpos = 0;
	} else if (res->outbufferpos > res->contentsize) {
		/* if we ended right at the end of the buffer/frame previous time */
		zstd_readheader(res);
		while (res->contentsize == 0) {
			if (feof(res->finput)) break;
			zstd_skipframe(res);
			zstd_readheader(res);
		}
		zstd_readframe(res);
		res->outbufferpos = 0;
	}
	read = 0;
	while (1) {
		char *match;
		writebuffer = res->outbuffer + res->outbufferpos;
		writesize = res->contentsize - res->outbufferpos;
		match = memchr(writebuffer, delim, writesize);
		if (match == NULL) {
			/* no delimiter found */
			if (ks_expand(str, writesize + 2) < 0) { return -3; }
			memcpy(str->s + str->l,writebuffer,writesize);
			str->l += writesize;
		} else {
			writesize = match - writebuffer;
			if (ks_expand(str, writesize + 2) < 0) { return -3; }
			memcpy(str->s + str->l,writebuffer,writesize);
			str->l += writesize;
			res->outbufferpos += writesize + 1;
			break;
		}
		zstd_readheader(res);
		while (res->contentsize == 0) {
			if (feof(res->finput)) break;
			zstd_skipframe(res);
			zstd_readheader(res);
		}
		if (feof(res->finput)) break;
		zstd_readframe(res);
	}
	read = str->l;
	if (read == 0 && feof(res->finput)) {return -1;}
	res->currentpos += read;
	if ( delim=='\n' && str->l>0 && str->s[str->l-1]=='\r' ) str->l--;
	str->s[str->l] = 0;
	return str->l <= INT_MAX ? (int) str->l : INT_MAX;
}

tbx_t *zst_tbx_index(ZSTDres *fp, int min_shift, const tbx_conf_t *conf)
{
    tbx_t *tbx;
    kstring_t str;
    int ret, first = 0, n_lvls, fmt;
    int64_t lineno = 0;
    uint64_t last_off = 0;
    tbx_intv_t intv;
    int64_t max_ref_len = 0;

    str.s = 0; str.l = str.m = 0;
    tbx = (tbx_t*)calloc(1, sizeof(tbx_t));
    if (!tbx) return NULL;
    tbx->conf = *conf;
    if (min_shift > 0) n_lvls = (TBX_MAX_SHIFT - min_shift + 2) / 3, fmt = HTS_FMT_CSI;
    else min_shift = 14, n_lvls = 5, fmt = HTS_FMT_TBI;
    while ((ret = zstd_tbx_getline(fp, '\n', &str)) >= 0) {
//fprintf(stderr,"-dbg- str=%*.*s\n",(int)str.l,(int)str.l,str.s);
        ++lineno;
        if (str.s[0] == tbx->conf.meta_char && fmt == HTS_FMT_CSI) {
            switch (tbx->conf.preset) {
                case TBX_SAM:
                    zst_adjust_max_ref_len_sam(str.s, &max_ref_len); break;
                case TBX_VCF:
                    zst_adjust_max_ref_len_vcf(str.s, &max_ref_len); break;
                default:
                    break;
            }
        }
        if (lineno <= tbx->conf.line_skip || str.s[0] == tbx->conf.meta_char) {
            last_off = zst_tbx_tell(fp);
            continue;
        }
//fprintf(stderr,"-dbg- first=%d\n",first);
        if (first == 0) {
            if (fmt == HTS_FMT_CSI) {
                if (!max_ref_len)
                    max_ref_len = (int64_t)100*1024*1024*1024; // 100G default
                n_lvls = zst_adjust_n_lvls(min_shift, n_lvls, max_ref_len);
            }
            tbx->idx = hts_idx_init(0, fmt, last_off, min_shift, n_lvls);
//fprintf(stderr,"-dbg- hts_idx_init\n");
            if (!tbx->idx) goto fail;
            first = 1;
        }
//fprintf(stderr,"-dbg- zst_get_intv str=%*.*s\n",(int)str.l,(int)str.l,str.s);
       ret = zst_get_intv(tbx, &str, &intv, 1);
//fprintf(stderr,"-dbg- ret=%d str=%*.*s\n",ret,(int)str.l,(int)str.l,str.s);
        if (ret < -1) goto fail;  // Out of memory
        if (ret < 0) continue; // Skip unparsable lines
//fprintf(stderr,"-dbg- intv.tid=%d intv.beg=%ld intv.end=%ld\n",intv.tid, intv.beg, intv.end);
        if (hts_idx_push(tbx->idx, intv.tid, intv.beg, intv.end,
                         zst_tbx_tell(fp), 1) < 0) {
            goto fail;
        }
    }
    if (ret < -1) goto fail;
    if ( !tbx->idx ) tbx->idx = hts_idx_init(0, fmt, last_off, min_shift, n_lvls);   // empty file
    if (!tbx->idx) goto fail;
    if ( !tbx->dict ) tbx->dict = kh_init(s2i);
    if (!tbx->dict) goto fail;
    if (hts_idx_finish(tbx->idx, zst_tbx_tell(fp)) != 0) goto fail;
    if (zst_tbx_set_meta(tbx) != 0) goto fail;
    free(str.s);
    return tbx;

 fail:
    free(str.s);
    tbx_destroy(tbx);
    return NULL;
}

int zst_tbx_index_build(const char *fn, int min_shift, const tbx_conf_t *conf)
{
    const char *fnidx = NULL;
    tbx_t *tbx;
    ZSTDres *fp;
    int ret;
//fprintf(stdout,"-dbg- build index\n");
    if ((fp = zstd_openfile((char *)fn)) == 0) return -1;
    // if ( n_threads ) bgzf_mt(fp, n_threads, 256);
    // if ( bgzf_compression(fp) != bgzf ) { bgzf_close(fp); return -2; }
    tbx = zst_tbx_index(fp, min_shift, conf);
    zstd_close(fp);
    if ( !tbx ) return -1;
    ret = hts_idx_save_as(tbx->idx, fn, fnidx, min_shift > 0? HTS_FMT_CSI : HTS_FMT_TBI);
    tbx_destroy(tbx);
    return ret;
}
