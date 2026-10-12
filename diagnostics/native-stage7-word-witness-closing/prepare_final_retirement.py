"""Adapt the reviewed finite whole-batch retirement to these three cases/stage."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
OLD=ROOT/'attestations/native-vulkan/24a8e14/native-rhs-layer-clock-execution/files'
records={}
for name in ('run_final_retirement.py','verify_final_retirement.ps1'):
    p=OLD/name;text=p.read_text();before={'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
    text=text.replace(", 'placement'",'').replace(",'placement'",'').replace("','placement'", "'")
    if name.endswith('.ps1'):
        text=text.replace('[string]$Bootstrap)', '[string]$Bootstrap,[Parameter(Mandatory=$true)][string]$NativeStage)')
        text=text.replace('finite-native-clock-control-retirement','finite-native-word-witness-retirement')
        text=text.replace('this whole WORK only','this whole WORK and unique native stage only')
        text=text.replace("-and $_.CommandLine.IndexOf($Stage,[StringComparison]::OrdinalIgnoreCase) -ge 0", "-and ($_.CommandLine.IndexOf($Stage,[StringComparison]::OrdinalIgnoreCase) -ge 0 -or $_.CommandLine.IndexOf($NativeStage,[StringComparison]::OrdinalIgnoreCase) -ge 0)")
        text=text.replace("Accessible consumer still refers to task WORK", "Accessible consumer still refers to task WORK or native stage")
    else:
        old="'-Bootstrap',win(cleanup_bootstrap)]"
        assert text.count(old)==1
        text=text.replace(old,"'-Bootstrap',win(cleanup_bootstrap),'-NativeStage',win(ROOT/read(HERE/'source.json')['stage_root'])]")
        text=text.replace("HERE/'verify_final_retirement.ps1',*[", "HERE/'verify_final_retirement.ps1',HERE/'source.json',*[")
    target=WORK/name;assert not target.exists();target.write_text(text)
    records[name]={'basis':str(p.relative_to(ROOT)),'basis_seal':before,'current_seal':{'bytes':target.stat().st_size,'sha256':hashlib.sha256(target.read_bytes()).hexdigest()}}
(WORK/'final-retirement-preparation.json').write_text(json.dumps({'records':records,'scope':'Original40s finite retirement and observer-stop recovery unchanged. Exactly3completedcase owner/child births and accessible wholeWORK/native-stage commandline references. No hidden/global/scheduled/handle absence claim.'},indent=2)+'\n');print(json.dumps(records,indent=2))
