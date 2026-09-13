// Private interpreter registers in workgroup memory. Each invocation owns one
// disjoint lane; there is no cross-lane dependency or barrier. The generated lane
// count keeps the largest program below Vulkan's 16 KiB shared-memory minimum.
groupshared uint retainedScratch[SIRIUS_RETAINED_REGISTERS * SIRIUS_RETAINED_TERMS * SIRIUS_RETAINED_LANES];
uint ScratchAddress(uint index, uint component, uint rowBase, uint prefix) {
    uint lane = (rowBase / (prefix + SIRIUS_RETAINED_TERMS * SIRIUS_RETAINED_REGISTERS)) % SIRIUS_RETAINED_LANES;
    return (index * SIRIUS_RETAINED_TERMS + component) * SIRIUS_RETAINED_LANES + lane;
}
uint ScratchRead(uint index, uint component, uint rowBase, uint prefix) {
    return retainedScratch[ScratchAddress(index, component, rowBase, prefix)];
}
void ScratchWrite(uint index, uint component, uint rowBase, uint prefix, uint value) {
    retainedScratch[ScratchAddress(index, component, rowBase, prefix)] = value;
}
