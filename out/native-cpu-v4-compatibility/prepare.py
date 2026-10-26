from pathlib import Path
import hashlib
basis=Path.cwd()/'out/native-cpu-local-identities-compatibility/prepare.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest()=='2f61022f4c9cdbb1493261ff5811b404d608ff63f5bb07492f07ab5ddf727c24'
source=basis.read_text()
assert source.count('bin/linux-gcc/native-cpu-local-identities-compatibility')==1
source=source.replace('bin/linux-gcc/native-cpu-local-identities-compatibility','bin/linux-gcc/native-cpu-v4-compatibility',1)
assert source.count('out/native-cpu-local-identities/provider/source-contract.json')==1
source=source.replace('out/native-cpu-local-identities/provider/source-contract.json','out/native-cpu-v4-provider/provider/source-contract.json',1)
assert source.count('Reuse exact genuine R fixture driver and renamed Q baseline object; both already compiled/verified. New AN raw3 with original integer identity stores before unchanged AM CPP calls plus unchanged U scalar3/pureF and combined original Iadapter/checked dense C++, one new link, zero CPP/old-math recompilations or objcopy operations. Original38 inputs/capacity40/187960words and scoped FP/ISA/admission checks unchanged. Historical R accounting failures stay historical.')==1
source=source.replace('Reuse exact genuine R fixture driver and renamed Q baseline object; both already compiled/verified. New AN raw3 with original integer identity stores before unchanged AM CPP calls plus unchanged U scalar3/pureF and combined original Iadapter/checked dense C++, one new link, zero CPP/old-math recompilations or objcopy operations. Original38 inputs/capacity40/187960words and scoped FP/ISA/admission checks unchanged. Historical R accounting failures stay historical.','Reuse exact genuine R fixture driver, renamed Q baseline and unchanged AN raw3/U scalar3/pureF. Only the8th merged adapter composes exact original v3 math with separately guarded same-eight-lane optional v4 math/baseline dispatch. One new link, zero C++/raw/host/SPIRV compilations. Original38 inputs/capacity40/187960words, counts/tolerances/FP/ISAeligibility/admission checks unchanged. Supporting finite compatibility only, no provider/performance/scientific/release admission.',1)
assert source.count('One fresh AN Transport compatibility workload, all38 unchanged cases/all187960words.')==1
source=source.replace('One fresh AN Transport compatibility workload, all38 unchanged cases/all187960words.','One fresh guarded-v4-composition Transport compatibility workload, all38 unchanged cases/all187960words.',1)
exec(compile(source,str(basis),"exec"))
