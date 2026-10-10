#ifdef SIRIUS_RETAINED_PAIRED_TRANSPORT
[[vk::binding(1, 0)]] RWStructuredBuffer<uint> outputs;
// Two independent row partitions share one workgroup. Five words remain
// explicit. The complete 459-register domain includes an output-backed spill.
groupshared uint retainedScratch[2 * 384 * 5];
groupshared uint retainedStatus[2];
uint RetainedPartition() {
    return spirv_asm { OpLoad $$uint result builtin(LocalInvocationIndex:uint); } / 32;
}
uint ScratchRead(uint index, uint component, uint rowBase, uint prefix) {
    if (index < 384)
        return retainedScratch[(RetainedPartition() * 384 + index) * 5 + component];
    return outputs[rowBase + prefix + index * 5 + component];
}
void ScratchWrite(uint index, uint component, uint rowBase, uint prefix, uint value) {
    if (index < 384)
        retainedScratch[(RetainedPartition() * 384 + index) * 5 + component] = value;
    else
        outputs[rowBase + prefix + index * 5 + component] = value;
}
void RejectRetainedProgram() {
    uint previous;
    InterlockedOr(retainedStatus[RetainedPartition()], 1, previous);
}
bool RetainedFailed() { return (retainedStatus[RetainedPartition()] & 1) != 0; }
bool AnyRetainedLive() {
    return (retainedStatus[0] & 0x80000001u) == 0x80000000u ||
           (retainedStatus[1] & 0x80000001u) == 0x80000000u;
}
uint RetainedLiveStatus(uint layers) { return 0x80000000u | (layers << 16); }
uint RetainedRounds() {
    uint a = (retainedStatus[0] & 0x80000001u) == 0x80000000u ?
             ((retainedStatus[0] >> 16) & 255) : 0;
    uint b = (retainedStatus[1] & 0x80000001u) == 0x80000000u ?
             ((retainedStatus[1] >> 16) & 255) : 0;
    return max(a, b);
}
#else
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
#endif
