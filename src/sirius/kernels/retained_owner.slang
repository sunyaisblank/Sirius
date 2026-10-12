// Optional portable high-row intermediates. Canonical references/maps are
// verified at provisioning and occupy an immutable prefix of this buffer.
import retained_owner_layout;
import retained_scalar;
import retained_triple;
#if !defined(SIRIUS_RETAINED_PORTABLE) || !defined(SIRIUS_RETAINED_NORMAL_SUM32)
#error Stage-seven ownership requires the portable normal-sum arithmetic route.
#endif
#if SIRIUS_RETAINED_EXECUTION_LANES != 16 || (!defined(SIRIUS_RETAINED_PARALLEL_TRANSPORT) && !defined(SIRIUS_RETAINED_PARALLEL_ENDPOINT))
#error Stage-seven ownership requires a parallel sixteen-invocation Transport or Endpoint.
#endif
[[vk::binding(2, 0)]] RWStructuredBuffer<uint> stageSevenOwner;
groupshared uint s7Enabled;
groupshared uint s7RowBase;
groupshared uint s7Family;
groupshared uint s7Invalid;
groupshared uint s7Captured;

// Control: version, capacity, generation low/high, arithmetic mode, grant words.
// Capacity-indexed grant bits follow control and shift every immutable offset.
// Row header: ready, reserved, family, generation low/high, reserved[3].
bool S7Shape(uint row, uint rows, out uint base) {
    base = 0;
    uint words, stride;
    stageSevenOwner.GetDimensions(words, stride);
    uint grants = (rows + 31u) / 32u;
    if (rows == 0 || row >= rows || words < S7PrefixWords + grants ||
        (words - S7PrefixWords - grants) % S7RowWords != 0 ||
        (words - S7PrefixWords - grants) / S7RowWords != rows) return false;
    if (stageSevenOwner[0] != 1 || stageSevenOwner[1] != rows ||
        stageSevenOwner[4] != 1 || stageSevenOwner[5] != grants ||
        (stageSevenOwner[2] | stageSevenOwner[3]) == 0) return false;
    base = S7PrefixWords + grants + row * S7RowWords;
    return true;
}
void S7ResetTransport(uint row, uint rows, uint lane) {
    if (lane != 0) return;
    s7Enabled = 0;
    s7RowBase = 0;
    s7Family = 0xffffffffu;
    s7Invalid = 0;
    s7Captured = 0;
    uint base;
    if (S7Shape(row, rows, base)) {
        s7Enabled = 1;
        s7RowBase = base;
        stageSevenOwner[base] = 0;
        stageSevenOwner[base + 1] = 0;
        stageSevenOwner[base + 3] = stageSevenOwner[2];
        stageSevenOwner[base + 4] = stageSevenOwner[3];
    }
}
void S7Store(uint address, RetainedTriple value) {
    stageSevenOwner[address] = RPBits(value.hi);
    stageSevenOwner[address + 1] = RPBits(value.lo);
    stageSevenOwner[address + 2] = RPBits(value.tail);
    stageSevenOwner[address + 3] = RPBits(value.error);
    stageSevenOwner[address + 4] = value.valid;
}
void S7Key(uint index, RetainedTriple value) {
    if (s7Enabled == 0 || (index >= 8 && index != 44)) return;
    uint slot = index == 44 ? 8 : index;
    S7Store(s7RowBase + S7KeyOffset + 5 * slot, value);
}
uint S7ValueCount() { return s7Family == 0 ? S7GeneralValues : S7SchwarzschildValues; }
uint S7Family(uint actualBase, uint rows, uint packetWords, uint generalWords) {
    uint first = 1 + rows * packetWords;
    return actualBase == first ? 0 :
           (actualBase == first + generalWords ? 1 : 0xffffffffu);
}
// Collective: every invocation in the row calls this after family selection.
// A miss changes only owner eligibility; it never touches retainedStatus.
void S7Recognize(StructuredBuffer<uint> actual, uint actualBase,
                 uint referenceBase, uint words, uint lane) {
    if (s7Enabled != 0)
        for (uint i = lane; i < words; i += SIRIUS_RETAINED_EXECUTION_LANES)
            if (actual[actualBase + i] != stageSevenOwner[referenceBase + stageSevenOwner[5] + i])
                InterlockedOr(s7Invalid, 1u);
    GroupMemoryBarrierWithGroupSync();
    if (lane == 0 && s7Invalid != 0) s7Enabled = 0;
    GroupMemoryBarrierWithGroupSync();
}
void S7TransportTable(StructuredBuffer<uint> actual, uint actualBase, uint rows, uint lane) {
    if (lane == 0) {
        s7Family = S7Family(actualBase, rows, 230, S7GeneralTransportTableWords);
        if (s7Family > 1) s7Enabled = 0;
    }
    GroupMemoryBarrierWithGroupSync();
    uint offset = s7Family == 0 ? S7GeneralTransportTableOffset : S7SchwarzschildTransportTableOffset;
    uint words = s7Family == 0 ? S7GeneralTransportTableWords : S7SchwarzschildTransportTableWords;
    S7Recognize(actual, actualBase, S7TransportReferenceOffset + offset, words, lane);
}
void S7Capture(uint ordinal, RetainedTriple value) {
    if (s7Enabled == 0) return;
    uint count = s7Family == 0 ? S7GeneralTransportInstructions : S7SchwarzschildTransportInstructions;
    if (ordinal >= count) return;
    uint map = s7Family == 0 ? S7GeneralTransportMapOffset : S7SchwarzschildTransportMapOffset;
    uint slot = stageSevenOwner[map + stageSevenOwner[5] + ordinal];
    if (slot == 0xffffffffu) return;
    if (slot >= S7ValueCount()) { InterlockedOr(s7Invalid, 1u); return; }
    S7Store(s7RowBase + S7ValuesOffset + 5 * slot, value);
    if (value.valid != 1) InterlockedOr(s7Invalid, 1u);
    InterlockedAdd(s7Captured, 1u);
}
// Called by lane zero only after the original final all-memory join and every
// original completion check/store. Invalid intermediates select generic replay.
void S7Publish() {
    if (s7Enabled == 0 || s7Invalid != 0 || s7Captured != S7ValueCount()) return;
    stageSevenOwner[s7RowBase + 2] = s7Family;
    stageSevenOwner[s7RowBase] = 1;
}
void S7EndpointTable(StructuredBuffer<uint> actual, RWStructuredBuffer<uint> packed,
                    uint actualBase, uint packedBase, uint row, uint rows, bool flat, uint lane) {
    if (lane == 0) {
        s7Enabled = 0;
        s7Invalid = 0;
        s7Captured = 0;
        s7RowBase = 0;
        s7Family = S7Family(actualBase, rows, 225, S7GeneralEndpointTableWords);
        uint base;
        if (!flat && s7Family < 2 && S7Shape(row, rows, base)) {
            s7RowBase = base;
            if (stageSevenOwner[base] == 1 &&
                (stageSevenOwner[S7ControlWords + row / 32u] & (1u << (row % 32u))) != 0 &&
                stageSevenOwner[base + 2] == s7Family &&
                stageSevenOwner[base + 3] == stageSevenOwner[2] &&
                stageSevenOwner[base + 4] == stageSevenOwner[3]) s7Enabled = 1;
        }
    }
    GroupMemoryBarrierWithGroupSync();
    if (s7Enabled != 0)
        for (uint key = lane; key < S7KeyCount; key += SIRIUS_RETAINED_EXECUTION_LANES) {
            uint index = key < 8 ? key : 44;
            for (uint word = 0; word < 5; ++word)
                if (packed[packedBase + 5 * index + word] !=
                    stageSevenOwner[s7RowBase + S7KeyOffset + 5 * key + word])
                    InterlockedOr(s7Invalid, 1u);
        }
    if (s7Enabled != 0)
        for (uint slot = lane; slot < S7ValueCount(); slot += SIRIUS_RETAINED_EXECUTION_LANES)
            if (stageSevenOwner[s7RowBase + S7ValuesOffset + 5 * slot + 4] != 1)
                InterlockedOr(s7Invalid, 1u);
    uint offset = s7Family == 0 ? S7GeneralEndpointTableOffset : S7SchwarzschildEndpointTableOffset;
    uint words = s7Family == 0 ? S7GeneralEndpointTableWords : S7SchwarzschildEndpointTableWords;
    S7Recognize(actual, actualBase, S7EndpointReferenceOffset + offset, words, lane);
}
bool S7Replay(uint ordinal, out RetainedTriple value) {
    value = RTInvalid();
    if (s7Enabled == 0) return false;
    uint count = s7Family == 0 ? S7GeneralEndpointInstructions : S7SchwarzschildEndpointInstructions;
    if (ordinal >= count) return false;
    uint map = s7Family == 0 ? S7GeneralEndpointMapOffset : S7SchwarzschildEndpointMapOffset;
    uint slot = stageSevenOwner[map + stageSevenOwner[5] + ordinal];
    if (slot >= S7ValueCount()) return false;
    uint base = s7RowBase + S7ValuesOffset + 5 * slot;
    RetainedTriple candidate;
    candidate.hi = RPFromBits(stageSevenOwner[base]);
    candidate.lo = RPFromBits(stageSevenOwner[base + 1]);
    candidate.tail = RPFromBits(stageSevenOwner[base + 2]);
    candidate.error = RPFromBits(stageSevenOwner[base + 3]);
    candidate.valid = stageSevenOwner[base + 4];
    if (candidate.valid != 1) return false;
    value = candidate;
    return true;
}
