// One ray owns a workgroup. A layer's disjoint expression results become
// visible at its barrier, before any input register can be reused.
// The generated register plan plus status fits the portable 16 KiB bound.
groupshared uint retainedScratch[SIRIUS_RETAINED_REGISTERS * SIRIUS_RETAINED_TERMS];
groupshared uint retainedStatus;
uint ScratchAddress(uint index, uint component, uint rowBase, uint prefix) {
    return index * SIRIUS_RETAINED_TERMS + component;
}
uint ScratchRead(uint index, uint component, uint rowBase, uint prefix) {
    return retainedScratch[ScratchAddress(index, component, rowBase, prefix)];
}
void ScratchWrite(uint index, uint component, uint rowBase, uint prefix, uint value) {
    retainedScratch[ScratchAddress(index, component, rowBase, prefix)] = value;
}
void RejectRetainedProgram() {
    uint previous;
    InterlockedOr(retainedStatus, 1, previous);
}
