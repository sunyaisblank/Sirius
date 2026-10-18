// One ray owns a workgroup. A layer's disjoint expression results become
// visible at its barrier, before any input register can be reused.
// The generated register plan plus status fits the portable 16 KiB bound.
#if SIRIUS_RETAINED_TERMS == 5 && !defined(SIRIUS_RETAINED_PORTABLE) && !defined(SIRIUS_RETAINED_FP64) && !defined(SIRIUS_RETAINED_STAGE7_OWNER)
groupshared uint retainedScratch[SIRIUS_RETAINED_REGISTERS * 4];
#else
groupshared uint retainedScratch[SIRIUS_RETAINED_REGISTERS * SIRIUS_RETAINED_TERMS];
#endif
groupshared uint retainedStatus;
#if SIRIUS_RETAINED_TERMS == 5 && !defined(SIRIUS_RETAINED_PORTABLE) && !defined(SIRIUS_RETAINED_FP64) && !defined(SIRIUS_RETAINED_STAGE7_OWNER)
uint ScratchAddress(uint index, uint component, uint rowBase, uint prefix) {
    return (component == 4u ? 3u : component) * SIRIUS_RETAINED_REGISTERS + index;
}
uint ScratchRead(uint index, uint component, uint rowBase, uint prefix) {
    uint word = retainedScratch[ScratchAddress(index, component, rowBase, prefix)];
    if (component == 3u) return word & 0x7fffffffu;
    if (component == 4u) return word >> 31u;
    return word;
}
void ScratchWrite(uint index, uint component, uint rowBase, uint prefix, uint value) {
    uint address = ScratchAddress(index, component, rowBase, prefix);
    if (component == 3u) {
        // The following same-invocation component-four write publishes validity.
        retainedScratch[address] = value & 0x7fffffffu;
    } else if (component == 4u) {
        retainedScratch[address] = (retainedScratch[address] & 0x7fffffffu) |
                                   ((value & 1u) << 31u);
    } else {
        retainedScratch[address] = value;
    }
}
#else
uint ScratchAddress(uint index, uint component, uint rowBase, uint prefix) {
    return index * SIRIUS_RETAINED_TERMS + component;
}
uint ScratchRead(uint index, uint component, uint rowBase, uint prefix) {
    return retainedScratch[ScratchAddress(index, component, rowBase, prefix)];
}
void ScratchWrite(uint index, uint component, uint rowBase, uint prefix, uint value) {
    retainedScratch[ScratchAddress(index, component, rowBase, prefix)] = value;
}
#endif
void RejectRetainedProgram() {
    uint previous;
    InterlockedOr(retainedStatus, 1, previous);
}
