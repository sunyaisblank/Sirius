from pathlib import Path
import hashlib,json,re,subprocess,sys
R=Path.cwd();W=Path(__file__).resolve().parent;B=R/'bin/linux-gcc/portable-checked-fp64';mode,folder=sys.argv[1:];folder=Path(folder)
def passed(name):
    owner=json.loads((W/name/'owner.json').read_text());assert owner['passed'] and owner['status']=='terminal' and not owner['remaining_owned_processes']
if mode=='build':
    for row in json.loads((W/'commands.json').read_text()):subprocess.run(row['argv'],cwd=row['cwd'],check=True)
    text=(B/'corpus.spvasm').read_text()
    assert sorted(re.findall(r'OpCapability (\S+)',text))==['Float64','Shader']
    assert set(re.findall(r'OpTypeInt (\d+) [01]',text))=={'32'}
    assert set(re.findall(r'OpTypeFloat (\d+)',text))=={'64'}
    assert 'OpExecutionMode %main LocalSize 64 1 1' in text
    assert 'OpFDiv' in text and ' Sqrt ' in text and 'OpConvertFToU' in text and 'OpUMulExtended' in text
    (folder/'compiler-scope.json').write_text(json.dumps({'capabilities':['Shader','Float64'],'integer_widths':[32],'float_widths':[64],'old_integer_fallbacks_retained':True,'no_new_float_execution_mode':True,'not_arithmetic_or_benefit_evidence':True},indent=2)+'\n')
elif mode=='host':
    passed('build')
    subprocess.run([str(B/'check_host'),str(W/'input.bin'),str(W/'expected.bin')],check=True)
else:
    assert mode=='run';passed('host')
    subprocess.run([str(B/'gpu_probe'),str(B/'corpus.spv'),str(W/'input.bin'),str(W/'expected.bin'),str(folder/'actual.bin'),'--dispatch'],check=True)
    assert (folder/'actual.bin').read_bytes()==(W/'expected.bin').read_bytes()
    (folder/'complete-output-agreement.json').write_text(json.dumps({'passed':True,'original_independent_records':174268,'complete_raw_words':522804,'mismatches':0,'actual_sha256':hashlib.file_digest((folder/'actual.bin').open('rb'),'sha256').hexdigest(),'scope':'Original independent scalar corpus; checked Float64 estimates with unchanged integer fallback. No integrated, benefit, frame or release qualification.'},indent=2)+'\n')
