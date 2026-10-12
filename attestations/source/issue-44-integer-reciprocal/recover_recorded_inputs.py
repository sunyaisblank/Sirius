"""Bind each recorded input version to exact archive bytes, Git or retained tools."""
from pathlib import Path
import hashlib,json,subprocess
r=Path.cwd();a=r/'attestations/software-vulkan/49bfb43/portable-integer-reciprocal';source='49bfb434b0a31df958572448acebf655eb192296';bundle=json.loads((a/'whole-byte-manifest.json').read_text());aliases={}
for name,v in bundle['files'].items():aliases.setdefault((v['bytes'],v['sha256']),[]).append('files/'+name)
tracked={p.decode() for p in subprocess.check_output(['git','ls-tree','-r','--name-only','-z',source]).split(b'\0') if p};records={};receipt_count=0
for p in (a/'files/out/portable-integer-reciprocal').rglob('inputs-before.json'):
 before=json.loads(p.read_text());after=json.loads(p.with_name('inputs-after.json').read_text());assert before==after,p;receipt_count+=1
 for name,v in before.items():records[(name,v['bytes'],v['sha256'])]=v
proof=[];external=a/'retained-input-snapshots';external.mkdir(exist_ok=True);counts={}
for (name,size,digest),v in sorted(records.items()):
 key=(size,digest);path=Path(name) if Path(name).is_absolute() else r/name;relative=str(path.relative_to(r)) if path.is_relative_to(r) else None
 if key in aliases:
  target=a/aliases[key][0];assert hashlib.sha256(target.read_bytes()).hexdigest()==digest;kind='archive';binding=aliases[key]
 elif relative in tracked:
  blob=subprocess.check_output(['git','rev-parse',source+':'+relative],text=True).strip();data=path.read_bytes() if path.is_file() and hashlib.sha256(path.read_bytes()).hexdigest()==digest else subprocess.check_output(['git','show',source+':'+relative]);assert len(data)==size and hashlib.sha256(data).hexdigest()==digest,(name,'git');kind='git';binding={'source_revision':source,'path':relative,'blob':blob}
 else:
  assert relative is None or not relative.startswith(('out/portable-integer-reciprocal/','bin/linux-gcc/portable-integer-reciprocal/')),(name,'owned input missing archive recovery')
  assert path.is_file(),name;data=path.read_bytes();assert len(data)==size and hashlib.sha256(data).hexdigest()==digest,(name,'retained mismatch');kind='retained';binding={'path':str(path),'resolved_path':str(path.resolve()),'bytes':size,'sha256':digest}
  # Preserve mutable installed external binaries in addition to their live tool identity.
  if relative is None:
   target=external/digest
   if not target.exists():target.write_bytes(data)
   assert hashlib.sha256(target.read_bytes()).hexdigest()==digest;binding['archive_copy']=str(target.relative_to(a))
 counts[kind]=counts.get(kind,0)+1;proof.append({'recorded_name':name,'bytes':size,'sha256':digest,'recovery_kind':kind,'binding':binding})
result={'passed':True,'source_revision':source,'recorded_input_versions':len(proof),'receipt_pairs':receipt_count,'recovery_counts':counts,'records':proof,'scope':'Every recorded admitted before/after input version recovered exactly. Unrecorded dynamic system dependencies and remote raw backup are not asserted.'};(a/'recorded-input-recovery.json').write_text(json.dumps(result,indent=2)+'\n');print({k:result[k] for k in ('passed','recorded_input_versions','receipt_pairs','recovery_counts')})
