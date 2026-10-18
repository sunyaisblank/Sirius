from pathlib import Path
import json,subprocess,sys
q=Path(__file__).resolve().parent
out=Path(sys.argv[1]);assert out.is_dir()
commands=json.loads((q/'commands.json').read_text());results=[]
for row in commands:
 argv=[x.replace('{folder}',str(out)) for x in row['argv']]
 result=subprocess.run(argv,capture_output=True,timeout=30)
 (out/(row['label']+'.stdout')).write_bytes(result.stdout)
 (out/(row['label']+'.stderr')).write_bytes(result.stderr)
 results.append(dict(label=row['label'],argv=argv,returncode=result.returncode))
 (out/'results.json').write_text(json.dumps(results,indent=2)+'\n')
 assert result.returncode==0,row['label']
print('Disassembled25changedpairs;50commands0;no compiler or provider')
