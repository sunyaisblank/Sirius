#include <dlfcn.h>
#include <errno.h>
#include <fcntl.h>
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static _Atomic unsigned sequence;
static _Atomic unsigned long long total_bytes;
static void output(int fd, const char *data, size_t size) {
    while (size) {
        ssize_t n = write(fd, data, size);
        if (n < 0 && errno == EINTR) continue;
        if (n <= 0) _exit(96);
        data += n; size -= (size_t)n;
    }
    if (close(fd)) _exit(96);
}

/* Temporary diagnostic: forward the actual pass unchanged, then serialize its
 * completed module. No production flags, IR attributes, or pass options change. */
void *LLVMRunPasses(void *module, const char *pipeline, void *target, void *options) {
    void *llvm = dlopen("/lib/x86_64-linux-gnu/libLLVM.so.20.1", RTLD_LAZY | RTLD_NOLOAD);
    if (!llvm) { fputs("optimizer observer: LLVM not already loaded\n", stderr); _exit(96); }
    void *(*run)(void *, const char *, void *, void *) = dlsym(llvm, "LLVMRunPasses");
    char *(*print)(void *) = dlsym(llvm, "LLVMPrintModuleToString");
    void (*dispose)(char *) = dlsym(llvm, "LLVMDisposeMessage");
    if (!run || !print || !dispose) { fputs("optimizer observer: missing LLVM API\n", stderr); _exit(96); }
    void *error = run(module, pipeline, target, options);
    if (error || !strstr(pipeline, "instcombine")) { dlclose(llvm); return error; }
    char *text = print(module);
    if (!text) _exit(96);
    size_t bytes = strnlen(text, 64u * 1024u * 1024u);
    unsigned id = atomic_fetch_add(&sequence, 1);
    if (bytes == 64u * 1024u * 1024u || id >= 32 ||
        atomic_fetch_add(&total_bytes, bytes) + bytes > 128u * 1024u * 1024u) _exit(96);
    const char *folder = getenv("SIRIUS_OPT_IR_FOLDER");
    char name[64], path[4096], record[128];
    if (!folder) _exit(96);
    snprintf(name, sizeof(name), "optimized-%02u.ll", id);
    if (snprintf(path, sizeof(path), "%s/%s", folder, name) >= (int)sizeof(path)) _exit(96);
    int fd = open(path, O_WRONLY | O_CREAT | O_EXCL, 0600);
    if (fd < 0) _exit(96);
    output(fd, text, bytes);
    dispose(text);
    if (snprintf(path, sizeof(path), "%s/optimizer-exports.tsv", folder) >= (int)sizeof(path)) _exit(96);
    int n = snprintf(record, sizeof(record), "%p\t%s\t%zu\n", module, name, bytes);
    if (n < 0 || n >= (int)sizeof(record)) _exit(96);
    fd = open(path, O_WRONLY | O_CREAT | O_APPEND, 0600);
    if (fd < 0) _exit(96);
    output(fd, record, (size_t)n);
    dlclose(llvm);
    return error;
}
