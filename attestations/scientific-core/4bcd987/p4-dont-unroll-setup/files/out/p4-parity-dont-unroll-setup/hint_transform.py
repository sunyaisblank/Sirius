"""Exact binary-only replacement of Unroll by DontUnroll; no tool execution."""
import collections
import hashlib
import struct

STOCK_SHA = "19991852c2bf8025aaa30b4cb3029e7bb980425d74de515efc4879df95eb0891"
STOCK_BYTES = 5562236


def derive(stock):
    if len(stock) != STOCK_BYTES or hashlib.sha256(stock).hexdigest() != STOCK_SHA:
        raise ValueError("exact full original ordinary module required")
    words = struct.unpack("<" + "I" * (len(stock) // 4), stock)
    if words[0] != 0x07230203 or words[1] != 0x00010500 or not words[3] or words[4]:
        raise ValueError("original SPIR-V header differs")
    changed = []
    opcodes = collections.Counter()
    masks = collections.Counter()
    output = bytearray(stock)
    at = 5
    while at < len(words):
        count, opcode = words[at] >> 16, words[at] & 0xffff
        if not count or at + count > len(words):
            raise ValueError("instruction framing exceeds original module")
        opcodes[opcode] += 1
        if opcode == 246:  # OpLoopMerge: merge ID, continue ID, loop control.
            if count < 4:
                raise ValueError("truncated OpLoopMerge")
            mask = words[at + 3]
            masks[mask] += 1
            if mask & 1:
                if mask != 1:
                    raise ValueError("original hint mask must be exactly one")
                after = 2
                struct.pack_into("<I", output, 4 * (at + 3), after)
                changed.append({"instruction_word": at, "word_index": at + 3,
                                "byte_offset": 4 * (at + 3), "opcode": 246,
                                "instruction_words": count, "before": mask,
                                "after": after, "xor": mask ^ after})
        at += count
    if at != len(words) or opcodes[246] != 5103 or len(changed) != 2494:
        raise ValueError("original loop/hint inventory differs")
    candidate = bytes(output)
    after_words = struct.unpack("<" + "I" * len(words), candidate)
    actual_diff = [i for i, pair in enumerate(zip(words, after_words)) if pair[0] != pair[1]]
    if actual_diff != [row["word_index"] for row in changed]:
        raise ValueError("a nonselected word changed")
    restored = bytearray(candidate)
    for row in changed:
        if row["xor"] != 3 or row["before"] != 1 or row["after"] != 2:
            raise ValueError("change exceeds the original mask-one to mask-two policy")
        struct.pack_into("<I", restored, row["byte_offset"], row["before"])
    if bytes(restored) != stock or after_words[:5] != words[:5]:
        raise ValueError("exact reversal/header preservation failed")
    proof = {"schema_version": 1, "kind": "exact-oploopmerge-unroll-to-dont-unroll",
             "stock": {"bytes": len(stock), "sha256": STOCK_SHA},
             "candidate": {"bytes": len(candidate), "sha256": hashlib.sha256(candidate).hexdigest()},
             "word_count": len(words), "instruction_count": sum(opcodes.values()),
             "opcode_histogram": {str(k): v for k, v in sorted(opcodes.items())},
             "loop_merge_count": opcodes[246], "original_mask_histogram": {str(k): v for k, v in sorted(masks.items())},
             "changed_word_count": len(changed),
             "changed_byte_count": sum(a != b for a, b in zip(stock, candidate)),
             "all_other_words_identical": True, "header_identical": True,
             "all_bits_except_unroll_and_dont_unroll_preserved": True, "dont_unroll_added": True,
             "only_original_mask_one_replaced_by_two": True,
             "exact_reversal_equals_stock": True,
             "bounds_ids_opcodes_modes_functions_precision_and_arithmetic_words_identical": True,
             "changes": changed}
    return candidate, proof
