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
