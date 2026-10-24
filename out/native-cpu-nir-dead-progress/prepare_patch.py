"""Emit the single NIR fixed-point cleanup; never modify Mesa or objects."""
from pathlib import Path
import difflib
import hashlib
import json

root = Path.cwd()
work = Path(__file__).resolve().parent
mesa = root / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
name = 'src/gallium/auxiliary/gallivm/lp_bld_nir.c'


def pin(data):
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def publish(path, data):
    if path.exists():
        assert path.read_bytes() == data, path
    else:
        path.write_bytes(data)


hashes = {
    name: '90947f4f085415ac463019affe1e2725e5dc83ec56e24e0b301cba57205337f5',
    'src/compiler/nir/nir_search.c': '862e9c41143674400690a5a5b4415dc86989a97ca661c1d4fbe56f15a436b9d3',
    'src/compiler/nir/nir_lower_packing.c': 'd00cba28403a900de6437c0d5b4d6d83ab54e5e3fa3cae6491d4621f1d5ba9c1',
    'src/compiler/nir/nir_opt_constant_folding.c': '1ae1f50f9b1e5946b3893d576fa5a43c1da2a1cfa747076ad7ec971b7944ad17',
    'src/compiler/nir/nir_opt_dce.c': '6724b5f4b2ccee2843446f877232fc18ccb991e639eff5206dfd2f3339f1461e',
}
inputs, sources = {}, {}
for relative, expected in hashes.items():
    p = mesa / relative
    data = p.read_bytes()
    assert pin(data)['sha256'] == expected, relative
    inputs[str(p.relative_to(root))] = pin(data)
    sources[relative] = data.decode()
old = sources[name]
anchor = '      NIR_PASS(progress, nir, nir_lower_subgroups, &subgroups_options);\n'
addition = (anchor +
    '      /* Rewrites can orphan pure SSA trees. Remove them before another\n'
    '       * fixed-point iteration optimizes/lower-packs those dead trees. */\n'
    '      NIR_PASS(progress, nir, nir_opt_dce);\n')
assert old.count(anchor) == 1
new = old.replace(anchor, addition, 1)
assert new.replace(addition, anchor, 1) == old
patch = ''.join(difflib.unified_diff(old.splitlines(True), new.splitlines(True),
                                    n=3, fromfile='a/' + name, tofile='b/' + name))

# These are implementation facts, not a claim to have recovered the live NIR.
search = sources['src/compiler/nir/nir_search.c']
assert 'if (instr->type == nir_instr_type_alu)\n            nir_instr_worklist_push_tail(worklist, instr);' in search
assert 'nir_def_rewrite_uses(&instr->def, ssa_val);' in search
assert 'nir_instr_remove(&instr->instr);\n   exec_list_push_tail(dead_instrs, &instr->instr.node);' in search
packing = sources['src/compiler/nir/nir_lower_packing.c']
assert 'nir_def_replace(&alu_instr->def, dest);\n\n   return true;' in packing
folding = sources['src/compiler/nir/nir_opt_constant_folding.c']
assert 'nir_def_replace(&alu->def, imm);\n   nir_instr_free(&alu->instr);\n\n   return true;' in folding
dce = sources['src/compiler/nir/nir_opt_dce.c']
assert 'case nir_instr_type_call:\n   case nir_instr_type_jump:\n      return true;' in dce
assert '!(info->flags & NIR_INTRINSIC_CAN_ELIMINATE)' in dce
assert 'mark_src_live(&nif->condition, defs_live);' in dce
assert 'nir_foreach_src(instr, mark_live_cb, defs_live);' in dce
assert 'nir_instr_remove(instr);\n         exec_list_push_tail(dead_instrs, &instr->node);' in dce

provider_contract = root / 'out/native-cpu-simd-endpoint-regions/provider/source-contract.json'
candidate = json.loads(provider_contract.read_text())
objects = candidate['native_objects']
assert len(objects) == 8
observation = root / 'out/native-cpu-simd-endpoint-regions/benefit/runtime/b1/user-pc-observation.json'
for p in (provider_contract, observation):
    inputs[str(p.relative_to(root))] = pin(p.read_bytes())
proof = {
    'changed_file': name, 'basis': pin(old.encode()), 'candidate': pin(new.encode()),
    'patch': pin(patch.encode()), 'source_inputs': inputs,
    'native_objects_reused_without_recompilation': objects,
    'new_native_objects': 0,
    'source_facts': [
        'The sampled constant-folding/algebraic/packing/texture/instruction-lowering pass combination occurs in lp_build_opt_nir main fixed-point loop.',
        'That loop currently has no dead-code cleanup; nir_lower_tex walks all instructions on every iteration even for texture-free shaders.',
        'nir_algebraic_impl enqueues every ALU without a liveness gate. Its replacement removes the old root, not its now-unused source cone.',
        'Constant folding and packing likewise replace ALU roots and report progress even if their replacement has no live consumer.',
        'One end-of-loop DCE transitively removes dead SSA cones so later iterations cannot continue rewriting those cones.',
        'nir_opt_dce marks all calls/jumps, non-eliminable intrinsics and branch conditions live, recursively marks their source definitions, then removes only unmarked instructions.',
        'DCE does not rewrite retained live FP instruction expressions or their float-control metadata. No native shader, field, RHS, trial, output store, whitelist, eligibility, guard, query count, compiler flag or iteration limit is changed.',
    ],
    'limits': [
        'The passive IP samples identify active passes, not the precise shader/module or live NIR expression cycle.',
        'The source proves possible dead-expression progress and a semantics-preserving cleanup; it does not prove this is the complete cause of the terminal timeout or guarantee convergence of any live rewrite cycle.',
        'Root must integrate the patch with a new private provider/cache identity and perform the original gates. Existing eight native objects remain byte-identical.',
    ],
    'patch_applied': False, 'compiler_build_or_shader_execution_launched': False,
}
publish(work / 'provider.patch', patch.encode())
publish(work / 'source-proof.json', (json.dumps(proof, indent=2) + '\n').encode())
assert (mesa / name).read_text() == old
print(json.dumps({'source_patch_written': True, 'changed_files': 1,
                  'new_native_objects': 0, 'patch_applied': False, 'compiled_or_executed': False}))
