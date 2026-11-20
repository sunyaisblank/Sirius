"""Prepare an ordinary stock-Vulkan renderer with the certified continuation policy."""

import hashlib
import json
from pathlib import Path
import shlex
import subprocess

ROOT = Path.cwd()
WORK = ROOT/'out/strict-backend-science'
BUILD = ROOT/'bin/linux-ci'
SOURCE = ROOT/'src/sirius/render/vulkan_renderer.cpp'
TARGET = 'src/sirius/app/sirius_render_tests'

def pin(path):
    data=path.read_bytes()
    return {'path':str(path.relative_to(ROOT)), 'bytes':len(data),
            'sha256':hashlib.sha256(data).hexdigest()}

original=SOURCE.read_text()
old='    auto spherical_escape = backend::RetainedTraceExecutor::SphericalEscape::Disabled;'
new='''    auto spherical_escape = device.Info().kind == backend::DeviceKind::kSoftware
                                ? backend::RetainedTraceExecutor::SphericalEscape::Certified
                                : backend::RetainedTraceExecutor::SphericalEscape::Disabled;'''
assert original.count(old)==1
proof=WORK/'certified-link'
proof.mkdir()
source=proof/'vulkan_renderer.cpp'
source.write_text(original.replace(old,new))
recipe,=[rule for rule in json.loads((BUILD/'compile_commands.json').read_text())
         if rule['file']==str(SOURCE)]
compile_argv=shlex.split(recipe['command'])
obj=proof/'vulkan_renderer.o'
compile_argv[compile_argv.index(str(SOURCE))]=str(source)
compile_argv[compile_argv.index('-o')+1]=str(obj)
subprocess.run(compile_argv,cwd=recipe['directory'],check=True)
selected=[]
for line in subprocess.check_output(['/usr/bin/ninja','-C',str(BUILD),'-t','commands',TARGET],
                                   text=True).splitlines():
    args=shlex.split(line)
    if '-o' in args and args[args.index('-o')+1]==TARGET:
        selected.append(args)
tokens,=selected
assert tokens[:2]==[':','&&']
link_argv=tokens[2:tokens.index('&&',2)]
exe=proof/'sirius_render_tests'
link_argv[link_argv.index('-o')+1]=str(exe)
link_argv.insert(link_argv.index('-o'),str(obj))
subprocess.run(link_argv,cwd=BUILD,check=True)
subprocess.run(['/usr/bin/cp','-al',str((BUILD/TARGET).parent/'resources'),
                str(proof/'resources')],check=True)
assert (proof/'resources/assets/Starfield.png').samefile(ROOT/'assets/Starfield.png')
record={'checkout':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        'purpose':'Exploratory ordinary provider comparison with one software route-policy change',
        'original_source':pin(SOURCE),'candidate_source':pin(source),'object':pin(obj),
        'executable':pin(exe),'compile_argv':compile_argv,'link_argv':link_argv}
(proof/'link.json').write_text(json.dumps(record,sort_keys=True,indent=2)+'\n')
