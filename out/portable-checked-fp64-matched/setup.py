from pathlib import Path
import hashlib,json,shutil,struct,zipfile
R=Path.cwd();W=Path(__file__).resolve().parent;B=R/'bin/linux-gcc/portable-checked-fp64-matched';(B/'candidate').mkdir(parents=True);(W/'headers').mkdir()
bundle=R/'attestations/qualification/0109ab46/int64-admission-and-divide64-rejected'
manifest=json.loads((bundle/'MANIFEST.json').read_text())
with zipfile.ZipFile(bundle/'diagnostics-and-trial.zip') as z:
    (W/'audited-observer.cpp').write_bytes(z.read('out/portable-divide64-matched/audited-observer.cpp'))
    build=z.read('out/portable-divide64-matched/build.py').decode().replace('portable-divide64-matched','portable-checked-fp64-matched').replace('"""Two matched observers, one new module, one replaced backend archive member."""','"""One reused exact observer, three unique modules, one replaced archive member."""').replace("'shared_original_observer_compiled_once':True","'shared_original_observer_reused_byte_exact':True").replace("'original_observer_compilations':1","'original_observer_compilations':0")
    (W/'build.py').write_text(build)
    prepare=z.read('out/portable-divide64-matched/prepare.py').decode().replace('portable-divide64-matched','portable-checked-fp64-matched')
    old="('/tests/backend/retained_compute_test.cpp',W/'audited-observer.cpp',B/'audited-observer.o',[])"
    assert old in prepare;prepare=prepare.replace(','+old,'')
    prepare=prepare.replace("R/'out/portable-divide64-admission'/n","R/'out/portable-checked-fp64'/n")
    prepare=prepare.replace("'bin/linux-gcc/portable-divide64-admission/candidate-endpoint-normal.spv'","'bin/linux-gcc/portable-checked-fp64-modules/checked-endpoint.spv','bin/linux-gcc/portable-checked-fp64-modules/checked-transport.spv','bin/linux-gcc/portable-checked-fp64-modules/checked-dopri_phase.spv'")
    prepare=prepare.replace('out/portable-divide64-jit/ordinary-divide-review.json','out/portable-checked-fp64-jit/ordinary-code-review.json').replace('out/portable-divide64-corpus/run/complete-output-agreement.json','out/portable-checked-fp64/run/complete-output-agreement.json')
    prepare=prepare.replace('Compile only one retained-compute candidate selector/object and one shared original observer; reuse all other current GCC objects. Relink original baseline and candidate, verify exactly one archive member changes, preserve all35 old arrays and add one unique optional64 Endpoint module. No Vulkan dispatch, integrated correctness, benefit or release claim.','Compile only one candidate three-stage selector/object; reuse the byte-exact original audited observer and all other GCC objects. Preserve all35 old arrays, add three uniquely named checked-FP64 modules, verify exactly one archive member changes. No Vulkan dispatch or integrated/benefit qualification.')
    prepare=prepare.replace("'new_cpp_units':2,'observer_compilations':1","'new_cpp_units':1,'observer_compilations':0")
    (W/'prepare.py').write_text(prepare)
member='bin/linux-gcc/portable-divide64-matched/audited-observer.o'
with zipfile.ZipFile(bundle/'actual-products.zip') as z:
    content=z.read(member);assert hashlib.sha256(content).hexdigest()==manifest['files'][member]['sha256'];(B/'audited-observer.o').write_bytes(content)
header=['#pragma once','#include "'+str(R/'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h')+'"','namespace sirius::backend::retained_program {']
modules=json.loads((R/'out/portable-checked-fp64-modules/compiler-result.json').read_text())['modules']
for kind,value in modules.items():
    data=(R/value['module']).read_bytes();assert hashlib.sha256(data).hexdigest()==value['sha256'];words=struct.unpack('<'+'I'*(len(data)//4),data)
    header.append('inline constexpr std::array<std::uint32_t, '+str(len(words))+'> k'+kind+'PortableNormalSumCheckedFp64Shader{{')
    for i in range(0,len(words),16):header.append(','.join(str(x)+'u' for x in words[i:i+16])+',')
    header.append('}};')
header.append('}');(W/'headers/retained_kernels.h').write_text('\n'.join(header)+'\n')
source=(R/'src/sirius/backend/retained_compute.cpp').read_text()
for kind in modules:
    old='normal_sum ? std::span(k'+kind+'PortableNormalSumShader)'
    assert source.count(old)==1
    source=source.replace(old,'normal_sum ? (device.Info().supports_fp64\n                         ? std::span<const std::uint32_t>(k'+kind+'PortableNormalSumCheckedFp64Shader)\n                         : std::span<const std::uint32_t>(k'+kind+'PortableNormalSumShader))')
(W/'retained_compute.cpp').write_text(source)
(W/'source-scope.json').write_text(json.dumps({'canonical_public_CXX_source':'0109ab462edade524667d20303507aa09adcbaf6','reused_actual_GCC_objects_revision':'442612671802d1ffb965c5a3e3a27f7897a78267','public_CXX_content_equal':True,'modules':modules,'selector_delta':'Only portable normal Endpoint/Transport/DopriPhase, in either existing product mode, use checked estimates when actual enabled Float64 is present. Otherwise original modules. All native/FMA/conservative/owner routes unchanged.','all_original35_arrays_and_tables_preserved':True,'no_ODR_replacement':True,'original_audited_observer_byte_exact':manifest['files'][member],'observer_compilations':0,'no_public_library_or_shader_modification':True},indent=2)+'\n')
