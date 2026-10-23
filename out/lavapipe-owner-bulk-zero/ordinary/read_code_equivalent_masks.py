"""Validate the saved Boolean-only optimizer deltas; reuse the complete code gate."""
from pathlib import Path
import hashlib
basis=Path(__file__).resolve().parent/'read_code.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest()=='5eb4d91076d7566965ce3734d593eb8c224e8f06cda3306ace3dfcaa371a3978'
source=basis.read_text()
addon=r'''
_original_normalized_clear = normalized_clear
_mask_proof_rows = []
def normalized_clear(ir, count, candidate):
    # No new input restriction. owner_partials is 0..7 for every defined i32
    # argument: values >=8 select0. The zero path skips the loop; otherwise
    # indices begin1..7 and terminate after7. The mask phi starts allones and
    # only inserts literal0 at a valid lane, so every mask bit is nonpoison.
    # Consequently select(C,maskBit,false) == C & maskBit, including poison C.
    # AND reassociation preserves the same operands/poison propagation. Only
    # these observed Boolean forms and one confined SSA bijection are admitted.
    mask_source = ''' + "'''" + '''  %owner_remaining = sub i32 %8, %5
  %owner_tail = icmp ult i32 %owner_remaining, 8
  %owner_partials = select i1 %owner_tail, i32 %owner_remaining, i32 0
  %.not = icmp eq i32 %owner_partials, 0''' + "'''" + r'''
    loop_source = ''' + "'''" + '''  %loop_counter2.0 = phi i32 [ %20, %loop_begin1 ], [ %owner_partials, %entry ]
  %mask.1 = phi <8 x i32> [ %19, %loop_begin1 ], [ splat (i32 -1), %entry ]
  %19 = insertelement <8 x i32> %mask.1, i32 0, i32 %loop_counter2.0
  %20 = add i32 %loop_counter2.0, 1
  %21 = icmp ugt i32 %20, 7
  br i1 %21, label %endif-block, label %loop_begin1''' + "'''" + r'''
    assert mask_source in ir and loop_source in ir
    assert '  br i1 %.not, label %endif-block, label %loop_begin1' in ir
    assert '  %mask.0 = phi <8 x i32> [ splat (i32 -1), %entry ], [ %19, %loop_begin1 ]' in ir
    replaced = 0
    if not candidate:
        if count == 3060:
            old = ''' + "'''" + '''  %43079 = and <8 x i1> %ssa_65680.not, %875
  %43080 = icmp ne <8 x i32> %mask.0, zeroinitializer
  %43081 = select <8 x i1> %43079, <8 x i1> %43080, <8 x i1> zeroinitializer''' + "'''" + r'''
            new = ''' + "'''" + '''  %43079 = icmp ne <8 x i32> %mask.0, zeroinitializer
  %43080 = and <8 x i1> %ssa_65680.not, %43079
  %43081 = and <8 x i1> %43080, %875''' + "'''" + r'''
            assert ir.count(old) == 1
            # Intermediate names have no escaping uses. Final43081's uses
            # remain exact under the complete remainder comparison below.
            rest = ir.replace(old,'')
            assert not re.search(r'%(?:43079|43080)(?![0-9])', rest)
        else:
            assert count == 3145
            old = ''' + "'''" + '''  %8817 = and <8 x i1> %reg329.0, %591
  %8818 = icmp ne <8 x i32> %mask.0, zeroinitializer
  %8819 = select <8 x i1> %8817, <8 x i1> %8818, <8 x i1> zeroinitializer
  %8820 = load i32, ptr %29, align 4
  %8821 = load ptr, ptr %28, align 8
  %8822 = ashr i32 %8820, 2
  %8823 = insertelement <8 x i32> poison, i32 %8822, i64 0
  %8824 = shufflevector <8 x i32> %8823, <8 x i32> poison, <8 x i32> zeroinitializer
  %8825 = zext nneg <8 x i32> %ssa_12917 to <8 x i64>
  %channel_ptr3721 = getelementptr i32, ptr %8821, <8 x i64> %8825
  %oob_cmp3722 = icmp ult <8 x i32> %ssa_12917, %8824
  %mask3723 = and <8 x i1> %oob_cmp3722, %8819''' + "'''" + r'''
            new = ''' + "'''" + '''  %8817 = icmp ne <8 x i32> %mask.0, zeroinitializer
  %8818 = and <8 x i1> %reg329.0, %8817
  %8819 = load i32, ptr %29, align 4
  %8820 = load ptr, ptr %28, align 8
  %8821 = ashr i32 %8819, 2
  %8822 = insertelement <8 x i32> poison, i32 %8821, i64 0
  %8823 = shufflevector <8 x i32> %8822, <8 x i32> poison, <8 x i32> zeroinitializer
  %8824 = zext nneg <8 x i32> %ssa_12917 to <8 x i64>
  %channel_ptr3721 = getelementptr i32, ptr %8820, <8 x i64> %8824
  %oob_cmp3722 = icmp ult <8 x i32> %ssa_12917, %8823
  %8825 = and <8 x i1> %8818, %oob_cmp3722
  %mask3723 = and <8 x i1> %8825, %591''' + "'''" + r'''
            assert ir.count(old) == 1
            rest = ir.replace(old,'')
            assert not re.search(r'%(?:881[7-9]|882[0-5])(?![0-9])', rest)
        ir = ir.replace(old,new)
        defs = dict(re.findall(r'^  (%[A-Za-z0-9_.]+) = (.*)$', ir, re.M))
        pattern = r'^  (%[A-Za-z0-9_.]+) = select <8 x i1> (%[A-Za-z0-9_.]+), <8 x i1> (%[A-Za-z0-9_.]+), <8 x i1> zeroinitializer$'
        def replace(match):
            nonlocal replaced
            if defs[match[3]] != 'icmp ne <8 x i32> %mask.0, zeroinitializer': return match[0]
            replaced += 1
            return f'  {match[1]} = and <8 x i1> {match[2]}, {match[3]}'
        ir = re.sub(pattern,replace,ir,flags=re.M)
        assert replaced == {3060:16,3145:27}[count], (count,replaced)
    _mask_proof_rows.append({'candidate':candidate,'vectors':count,'simple_selects_canonicalized':replaced,
        'mask_full_definition_verified':True,'exact_local_reassociation_blocks_verified':not candidate,
        'scope':'Source/IR equivalence reasoning, not machine-checked proof or numerical acceptance.'})
    return _original_normalized_clear(ir,count,candidate)
'''
# Addon is deliberately limited to these pinned saved IR variants.
assert source.count('\ncompiler=')==1
source=source.replace('\ncompiler=', '\n'+addon+'\ncompiler=')
source=source.replace("'structural_premise_observed':True,'rows':rows,", "'structural_premise_observed':True,'rows':rows,'Boolean_mask_equivalence_review':_mask_proof_rows,")
source=source.replace('Remaining identifiers, operations, constants, flags, types and control flow exact.', 'After separately declared verified Boolean-mask identities and one confined SSA bijection, remaining identifiers, operations, constants, flags, types and control flow exact.')
exec(compile(source,str(basis),'exec'))
