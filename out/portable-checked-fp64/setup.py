from pathlib import Path
import json

R=Path.cwd();W=Path(__file__).resolve().parent;B=R/'bin/linux-gcc/portable-checked-fp64';B.mkdir()
rows=json.loads((W/'commands.json').read_text())
for row in rows:
    row['argv']=[x.replace('portable-divide64-corpus','portable-checked-fp64') for x in row['argv'] if x!='-DSIRIUS_PB32_DIVIDE64=1']
    row['argv']=[str(W) if 'portable-divide64-admission/candidate/src/sirius/kernels' in x else x for x in row['argv']]
rows.insert(0,{'cwd':str(R),'argv':['/usr/bin/g++-14','-std=c++2c','-O2','-g','-Wall','-Wextra','-Wpedantic','-Werror','-fsanitize=undefined,float-cast-overflow','-fno-sanitize-recover=all',str(W/'check_host.cpp'),'-o',str(B/'check_host')]})
(W/'commands.json').write_text(json.dumps(rows,indent=2)+'\n')
p=W/'gpu_probe.cpp';source=p.read_text().replace('!info.supports_int64','!info.supports_fp64').replace('actual supported Int64 software factory required','actual supported/enabled Float64 software factory required');p.write_text(source)
