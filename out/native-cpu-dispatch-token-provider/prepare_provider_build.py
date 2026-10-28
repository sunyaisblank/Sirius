"""Gate the same dry-run inventory and original bounded build ownership."""
from pathlib import Path
from authority import source, replace
from component_gate import component, contract_inputs, doc, terminal

source('out/native-cpu-inline-private-provider/prepare_provider_build.py', '1de3290baa4a3ab1d083b0340399baf73520613298e6e468a480bc3f6f8653b0')

source('out/native-cpu-packed-math-provider/prepare_provider_build.py', '016702edec2097b9b5495798a831c62470f6482bd557b7755540ce4a03b1241a')

source('out/native-cpu-v4-raw-provider/prepare_provider_build.py', '9504a012ea8ddc339c1d1dbf3f48a75ba9ea458f15f334623a11a6a684fa7ed3')

component()
P = Path(__file__).resolve().parent / 'provider'
contract_inputs(doc(P / 'source-contract.json'))
terminal(P / 'configure/owner.json')
source('out/native-cpu-v4-provider/prepare_provider_build.py', 'ccfd67ef343bf6c24fc141ae16b6601438c34c52063469cbef766ba59bd5718f')
source('out/native-cpu-endpoint-header-provider/prepare_provider_build.py', '5ebd23db67537967c89d562e22c43e89eb586a27411df385d1d35acbda34d4c2')
body = source('out/native-cpu-slang-bridge/prepare_provider_build.py', '5a91184fb45804da505f76552448fdb617a8eff6220909c5351bb031eee531d1')
body = replace(body, 'linking3 sealed same-source native CPU objects',
    'linking8 native rows with seven sealed inputs in the eighth, new guarded/admitted scalar3 and two exact header branches')
# The predecessor's exact planned outputs plus the single newly dependent C
# object are the complete allowed inventory; unknown commands stop before build.
plan_authority = Path.cwd() / 'out/native-cpu-inline-private-provider/provider/incremental-plan.json'
import hashlib, json
assert {'bytes': plan_authority.stat().st_size, 'sha256': hashlib.sha256(plan_authority.read_bytes()).hexdigest()} == {'bytes': 215726, 'sha256': '2ee9fe0b0ce2abe3777b90fb914868141c07eeeac1240d4ea814b706a3523311'}
expected_dirty_outputs = set(['bin/linux-gcc/lavapipe-acqrel-ordering/src/gallium/drivers/llvmpipe/libllvmpipe.a', 'bin/linux-gcc/lavapipe-acqrel-ordering/src/gallium/drivers/llvmpipe/libllvmpipe.a.p/lp_screen.c.o', 'bin/linux-gcc/lavapipe-acqrel-ordering/src/gallium/drivers/llvmpipe/libllvmpipe.a.p/lp_texture.c.o', 'bin/linux-gcc/lavapipe-acqrel-ordering/src/gallium/frontends/lavapipe/liblavapipe_st.a', 'bin/linux-gcc/lavapipe-acqrel-ordering/src/gallium/frontends/lavapipe/liblavapipe_st.a.p/lvp_device.c.o', 'bin/linux-gcc/lavapipe-acqrel-ordering/src/gallium/targets/lavapipe/libvulkan_lvp.so', 'bin/linux-gcc/lavapipe-acqrel-ordering/src/git_sha1.h', 'bin/linux-gcc/lavapipe-acqrel-ordering/src/loader/libloader.a', 'bin/linux-gcc/lavapipe-acqrel-ordering/src/loader/libloader.a.p/loader.c.o', 'bin/linux-gcc/lavapipe-acqrel-ordering/src/gallium/drivers/llvmpipe/libllvmpipe.a.p/lp_state_cs.c.o'])
body = replace(body, "target = B / 'src/gallium/targets/lavapipe/libvulkan_lvp.so'",
    "assert mutable == expected_dirty_outputs, ('Review changed dirty-output inventory before build', sorted(mutable ^ expected_dirty_outputs))\n" +
    "assert len(commands) == len(expected_dirty_outputs)\n" +
    "target = B / 'src/gallium/targets/lavapipe/libvulkan_lvp.so'")
exec(compile(body, str(Path(__file__)), 'exec'))
