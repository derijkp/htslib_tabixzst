struct ZSTF {
    // Reserved bits should be written as 0; read as "don't care"
    unsigned errcode:16, reserved:1, is_write:1, no_eof_block:1, is_be:1;
    signed compress_level:9;
    unsigned last_block_eof:1, is_compressed:1, is_gzip:1;
    int cache_size;
    int block_length, block_clength, block_offset;
    int64_t block_address, uncompressed_address;
    void *uncompressed_block, *compressed_block;
    bgzf_cache_t *cache;
    struct hFILE *fp; // actual file handle
    struct bgzf_mtaux_t *mt; // only used for multi-threading
    bgzidx_t *idx;      // BGZF index
    int idx_build_otf;  // build index on the fly, set by bgzf_index_build_init()
    struct z_stream_s *gz_stream; // for gzip-compressed files
    int64_t seeked;     // virtual offset of last seek
};
#ifndef HTS_BGZF_TYPEDEF
typedef struct BGZF BGZF;
#define HTS_BGZF_TYPEDEF
#endif

int zst_tbx_index_build(const char *fn, int min_shift, const tbx_conf_t *conf);
