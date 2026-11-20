"""Derive the SIMD16 Transport bounds checks from the sealed SIMD8 proof."""

from pathlib import Path

from resume_raw import R, pin


BASE = R / 'out/native-cpu-strict-ir-linking/closed_transport_context.cpp'
BASE_SHA = '7ac30679d18583c81e85e4aacd66f21d67c776ffaf9c659be936fa2fa8537677'


def derive(output):
    if pin(BASE)['sha256'] != BASE_SHA:
        raise RuntimeError('The reviewed Transport context proof changed')
    source = BASE.read_text()
    start = source.index('constexpr StringLiteral DomainKey =')
    owned = source.index('AllocaInst *ownedAllocation(', start)
    reset = source.index('unsigned expandReset(', owned)
    end = source.index('struct ClosedContext :', reset)
    body = source[start:owned] + source[reset:end]
    edits = (
        ('73472 + 32 * (3 * K + J / 2) : 160 + 160 * (3 * K + J / 2)',
         '146944 + 64 * (3 * K + J / 2) : 320 + 320 * (3 * K + J / 2)', 1),
        ('73472 + 32 * K : 160 + 160 * K', '146944 + 64 * K : 320 + 320 * K', 1),
        ('Offset->ugt(3968)', 'Offset->ugt(7936)', 2),
        ('Offset->getZExtValue() % 32', 'Offset->getZExtValue() % ColumnBytes', 2),
        ('initialized4000B field', 'initialized8000B field', 2),
        ('== ContextBytes', '== Transport.ContextBytes', 1),
        ('O >= 88192 && O < 92192', 'O >= 176384 && O < Transport.ContextBytes', 1),
        ('Entry.getInt32Ty(), 8)', 'Entry.getInt32Ty(), Lanes)', 1),
        ('"closed.coefficient.private4000"', '"wide.coefficient.private8000"', 1),
        ('Init.getInt8(0), 4000, Align(64)', 'Init.getInt8(0), 8000, Align(64)', 1),
        ('O - 88192', 'O - 176384', 1),
        ('Mask->getNumElements() != 8', 'Mask->getNumElements() != Lanes', 1),
        ('isolated_coefficient_bytes=4000', 'isolated_coefficient_bytes=8000', 1),
    )
    original = body
    for before, after, count in edits:
        if body.count(before) != count or after in body:
            raise RuntimeError(f'Transport proof substitution no longer matches {before}')
        body = body.replace(before, after)
    inverse = body
    for before, after, _ in reversed(edits):
        inverse = inverse.replace(after, before)
    if inverse != original:
        raise RuntimeError('Transport proof no longer inverts byte-for-byte')
    with Path(output).open('x') as stream:
        stream.write(body)
    return pin(output)
