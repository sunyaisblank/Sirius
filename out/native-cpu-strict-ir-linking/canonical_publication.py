"""Remove only escaping canonical publication records, preserving word order."""
from pathlib import Path
import hashlib,json,re,sys

R=Path(__file__).resolve().parents[2];W=Path(__file__).resolve().parent
B=R/'bin/linux-gcc/native-cpu-canonical-metadata/ispc-source'
STAGES=dict(dopri_phase=200,transport=80,endpoint=240)
FIELDS=['hi_0','lo_0','tail_0','error_0','valid_0']


def pin(data):return dict(bytes=len(data),sha256=hashlib.sha256(data).hexdigest())


def derive(stage,source):
 before='#define sirius_simd_'+stage+'_packet sirius_simd_canonical_'+stage+'_packet'
 after=before.replace('sirius_simd_canonical_','sirius_simd_publication_')
 assert source.count(before)==1 and after not in source
 # Original conversions are pure representation casts, including NaN payloads.
 assert re.search(r'static float RPFromBits_0\(uint32_t bits_0\)\s*\{\s*#line [^\n]+\n\s*return \(U32_asfloat\(\(bits_0\)\)\);\s*\}',source)
 assert re.search(r'static uint32_t RPBits_0\(float value_0\)\s*\{\s*#line [^\n]+\n\s*return \(F32_asuint\(\(value_0\)\)\);\s*\}',source)
 assert 'static uint32_t F32_asuint(float x) { return intbits(x); }' in source
 assert 'static float U32_asfloat(uint32_t x) { return floatbits(x); }' in source
 begin=source.index('    /* BEGIN pinned canonical metadata specialization. */\n')
 end=source.index('    /* END pinned canonical metadata specialization; exact fallback follows. */\n')
 fast=source[begin:end];edits=[]
 pattern=re.compile(r'(?P<pad> +)\{\n(?P=pad)    RetainedTriple_0 SiriusRoot;\n(?P<loads>(?:(?P=pad)    SiriusRoot\.[^\n]+\n){5})(?P=pad)    WriteStored_0\((?P<address>[^,\n]+), &SiriusRoot, (?P<context>\w+)\);\n(?P=pad)\}')
 def replace(m):
  context=m['context'];words=[]
  for field,line in zip(FIELDS,m['loads'].splitlines()):
   rhs=line.strip().removeprefix('SiriusRoot.'+field+' = ').removesuffix(';')
   if field!='valid_0':
    assert rhs.startswith('RPFromBits_0(') and rhs.endswith(')');rhs=rhs[len('RPFromBits_0('):-1]
   matched=re.fullmatch(re.escape(context)+r'->retainedScratch_0\[(\d+)U\]',rhs);assert matched,line
   words.append(int(matched[1]))
  assert words==list(range(words[0],words[0]+5)) and words[0]%5==0
  pad=m['pad'];lines=[pad+'{']
  for i,word in enumerate(words):lines.append(pad+f'    const uint32_t SiriusPublishedWord{i} = {context}->retainedScratch_0[{word}U];')
  for i in range(5):
   address=m['address'] if i==0 else m['address']+f' + {i}U'
   lines.append(pad+f'    *(&OUTPUT_WORD({context}, {address})) = SiriusPublishedWord{i};')
  lines.append(pad+'}');new='\n'.join(lines);edits.append((m[0],new));return new
 changed=pattern.sub(replace,fast)
 assert len(edits)==STAGES[stage],(stage,len(edits))
 assert 'WriteStored_0(' not in changed
 candidate=(source[:begin]+changed+source[end:]).replace(before,after,1)
 inverse=candidate.replace(after,before,1)
 for old,new in reversed(edits):
  assert new in inverse
  inverse=inverse.replace(new,old,1)
 assert inverse==source
 assert candidate.count('EvaluateOneStatic_0<')==source.count('EvaluateOneStatic_0<')
 return candidate,dict(publication_records_removed=len(edits),all5_source_words_loaded_before_any_store=True,
  same_original_five_OUTPUT_WORD_lvalues_in_order=True,bit_pattern_identity=True,complete_source_inverse_exact=True,
  original_fallback_unchanged=True,compiler_or_speed_claim=False)


def main():
 assert sys.argv[1:]==['derive']
 basis=json.loads((W/'canonical_metadata_contract.json').read_text())
 records={};generated={}
 for stage in STAGES:
  original=(B/(stage+'.ispc')).read_bytes();assert pin(original)==basis['derived_sources'][stage]
  candidate,proof=derive(stage,original.decode());target=B/(stage+'.publication.ispc')
  if target.exists():assert target.read_text()==candidate
  records[stage]=dict(input=pin(original),output=dict(path=str(target.relative_to(R)),**pin(candidate.encode())),source_proof=proof)
  generated[target]=candidate
 for target,candidate in generated.items():
  if not target.exists():
   with target.open('x') as s:s.write(candidate)
 destination=R/'out/native-cpu-canonical-metadata/publication-source.json';assert not destination.exists()
 result=dict(generator=dict(path=str(Path(__file__).relative_to(R)),**pin(Path(__file__).read_bytes())),stages=records,
  precondition='Existing same-thread positive immutable-program admission; internal scratch is disjoint from global output. All required roots validated before first publication. No source buffer/layout/allocation/count/arithmetic/rounding/tolerance change.',
  scope='Only inline the literal5word representation publication of the already-generated canonical fast branch. No compilation, runtime, numerical or speed qualification.')
 destination.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
