"""Preserve the pinned public sources used for the RDP detailed-target boundary."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,urllib.request
ROOT=Path.cwd();W=Path(__file__).resolve().parent
PAL='c5e800072a32f68b6ccc4422936d96167c6e0728';RDP='d63073a854c8ccc4b10edbf27559bb698c5685b3'
sources=[
 ('GPUOpen-Drivers/pal',PAL,'src/gpuUtil/gpuPerfExperimentTraceSource.h',6543,'ed219731c2c45ab994e15bdf2656ff3d6446677d17fe6a992f2925dafdf64466'),
 ('GPUOpen-Drivers/pal',PAL,'src/gpuUtil/gpuPerfExperimentTraceSource.cpp',27415,'74d648fbc7a09d4090688c435aa36449c837dff3e3be3cb03dc3380804b3244d'),
 ('GPUOpen-Drivers/pal',PAL,'src/gpuUtil/gpaSession.cpp',214614,'b1724567a4c8df797805a24b8954e194d86958f992b0a5ac1005575ae3fd65e2'),
 ('GPUOpen-Drivers/pal',PAL,'src/core/hw/gfxip/gfx9/gfx9PerfExperiment.cpp',203956,'f35290add15042c55ce2edfa8cc33af2e94ee4cdd284d9ed0739aec664ba960b'),
 ('GPUOpen-Drivers/pal',PAL,'inc/core/palPerfExperiment.h',30389,'209c06acad87872a838c6363b6aa7ce410d29eb25aca4eb9d14a9eed5ceabff9'),
 ('GPUOpen-Drivers/pal',PAL,'src/core/settings_core.json',66130,'b1f6278ccb6af76616acbc3be960bb81339a10dc196195bde48083c70bfc9437'),
 ('GPUOpen-Drivers/pal',PAL,'src/gpuUtil/asicInfoTraceSource.h',9747,'8513e70be0b6e053f44680400da6f485a8de8a0de65756e9bf4c675ef540d5aa'),
 ('GPUOpen-Drivers/pal',PAL,'src/gpuUtil/asicInfoTraceSource.cpp',12279,'c8b198e04d1036542f9b45ef246dfa3d3babc050dd4fc4bfb87648e3cb1c98d5'),
 ('GPUOpen-Tools/radeon_developer_panel',RDP,'source/api/RdpCaptureApi.h',57173,'20cd77a60584ae9f75aa47ace0fda4787918828f6cc566c1a2ee162175d71e52'),
 ('GPUOpen-Tools/radeon_developer_panel',RDP,'source/api/capture_params.cpp',4362,'c6e74eae1cf0a723d6634ebbeb691cfcca9a01b7af5022eda634afc46fdcc745'),
 ('GPUOpen-Tools/radeon_developer_panel',RDP,'source/trace/inc/rgp_trace_source.h',4993,'b357a69eb656a30c53cde51cd0aa67cf40a1aaa0afe128b2e0a78ec0d01f1d65'),
 ('GPUOpen-Tools/radeon_developer_panel',RDP,'source/trace/src/rgp/rgp_ubertrace_client.cpp',15668,'a0b03470e265a57f1809b6e74ec3ab1c991237484bd79dfbce289a347f0881cc'),
 ('GPUOpen-Tools/radeon_developer_panel',RDP,'source/trace/src/base_trace_source/ubertrace_params.h',15679,'b5a17af00d1509033cd88528797705ba37da3d3be9174ca759e1c99a0ab8011c'),
 ('GPUOpen-Tools/radeon_developer_panel',RDP,'source/cli/capture_cli.cpp',137180,'97cd53f709f425a65ab8e8ecce99eba6c47addf2b2475a0892046ecd9e22bddf')
]
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
records=[]
for repo,rev,source,size,digest in sources:
 url='https://raw.githubusercontent.com/'+repo+'/'+rev+'/'+source
 with urllib.request.urlopen(url,timeout=20) as response:data=response.read(1024**2)
 assert len(data)==size and hashlib.sha256(data).hexdigest()==digest
 p=W/('lane-'+('pal-' if repo.endswith('/pal') else 'rdp-')+Path(source).name);assert not p.exists();p.write_bytes(data)
 records.append({'path':str(p.relative_to(ROOT)),'url':url,'revision':rev,'repository_path':source,**seal(p)})
p=W/'lane-source-preparation.json';assert not p.exists();p.write_text(json.dumps({'files':records,'input_seals':{str(Path(__file__).relative_to(ROOT)):seal(Path(__file__))},'scope':'Pinned public sources used to settle supported detailed-target controls; no native/device execution or mutation.'},indent=2)+'\n')
print(json.dumps({'files':len(records),'bytes':sum(r['bytes'] for r in records)}))
