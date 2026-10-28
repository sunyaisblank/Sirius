"""Retain the original bounded dirty-output inventory; headers are unchanged."""
from pathlib import Path
from authority import source, replace
from component_gate import component, contract_inputs, doc, pin, terminal

component()
P = Path(__file__).resolve().parent / 'provider'
contract_inputs(doc(P / 'source-contract.json'))
terminal(P / 'configure/owner.json')
plan_authority = Path.cwd() / 'out/native-cpu-dispatch-token-provider/provider/incremental-plan.json'
assert pin(plan_authority) == dict(bytes=219780,
    sha256='49e3068ced436509fa59ebe9ed9f23c66d383160c6ea297e2d78132254ee9ae5')
expected_dirty_outputs = set(doc(plan_authority)['mutable_outputs'])
expected_dirty_outputs.remove('bin/linux-gcc/lavapipe-acqrel-ordering/src/gallium/drivers/llvmpipe/libllvmpipe.a.p/lp_state_cs.c.o')
assert len(expected_dirty_outputs) == 9
body = source('out/native-cpu-slang-bridge/prepare_provider_build.py',
    '5a91184fb45804da505f76552448fdb617a8eff6220909c5351bb031eee531d1')
body = replace(body, 'linking3 sealed same-source native CPU objects',
    'linking8 native rows with an eleven-input shared-value eighth object; both AT headers/first7 rows unchanged')
body = replace(body, "target = B / 'src/gallium/targets/lavapipe/libvulkan_lvp.so'",
    "assert mutable == expected_dirty_outputs, ('Review changed dirty-output inventory before build', sorted(mutable ^ expected_dirty_outputs))\n" +
    "assert len(commands) == len(expected_dirty_outputs)\n" +
    "target = B / 'src/gallium/targets/lavapipe/libvulkan_lvp.so'")
exec(compile(body, str(Path(__file__)), 'exec'))
