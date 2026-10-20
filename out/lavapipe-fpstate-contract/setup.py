from pathlib import Path
import json,shlex,shutil,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;B=R/'bin/linux-gcc/lavapipe-fpstate-contract';B.mkdir()
(W/'word_io.h').write_bytes(subprocess.check_output(['git','show','c41bc441:out/portable-checked-fp64/word_io.h']))
for name in ['run_owned.py','owned_guard.py','input.bin','expected.bin','portable_binary32.h','oracle-binding.json','corpus.slang','capture_jit.py','capture_features.py','read_features.py','feature-layout.json','properties_capture.py','retained-endpoint-arguments.json']:
    shutil.copyfile(R/'out/lavapipe-native-fp32'/name,W/name)
source=(R/'out/lavapipe-native-fp32/gpu_probe.cpp').read_text()
source=source.replace('#include <array>','#include <array>\n#include <xmmintrin.h>')
needle='int main(int argc, char** argv) {';assert source.count(needle)==1
source=source.replace(needle,'''// A Vulkan client may use a different host rounding/denormal environment.
// The shader explicitly requests the advertised RTE32 mode. Keep the caller
// state private to this process and restore it after all device resources.
struct CallerFpState {
    unsigned saved = _mm_getcsr();
    CallerFpState() {
        const unsigned requested = (saved & ~0x6000u) | 0x4000u | 0x8040u;
        _mm_setcsr(requested);
        std::cout << std::unitbuf << "{\\"caller_fpstate_before\\":" << saved
                  << ",\\"caller_fpstate_requested\\":" << requested
                  << ",\\"caller_fpstate_observed\\":" << _mm_getcsr() << "}\\n";
    }
    ~CallerFpState() { _mm_setcsr(saved); }
};
int main(int argc, char** argv) {
    CallerFpState caller;
''')
(W/'gpu_probe.cpp').write_text(source)
rows=json.loads((R/'out/lavapipe-native-fp32/original-probe-commands.json').read_text())[:3]
rows=[row for row in rows if '/usr/bin/g++-14'==row['argv'][0] and '-fsanitize=undefined,float-cast-overflow' not in row['argv']]
assert len(rows)==2
for row in rows:
    row['argv']=[arg.replace('portable-checked-fp64','lavapipe-fpstate-contract') for arg in row['argv']]
(W/'commands.json').write_text(json.dumps(rows,indent=2)+'\n')
for name in ['native.spv','native.spvasm','float-layout.json']:
    shutil.copyfile(R/'bin/linux-gcc/lavapipe-native-fp32'/name,B/name)
s=(R/'out/lavapipe-native-fp32/run_capture.py').read_text().replace('lavapipe-native-fp32','lavapipe-fpstate-contract')
s=s.replace('Default denormal behavior','Default denormal behavior under a non-default legal host FP environment')
(W/'run_capture.py').write_text(s)
