from pathlib import Path
import hashlib,json,re
work=Path(__file__).resolve().parent
source=Path('tests/backend/retained_compute_test.cpp').read_text()
cls=source[source.index('class PreparationProbeDevice final'):source.index('\nvoid ExpectPhysicalStatsEqual')]
tests=source[source.index('void ModelAcceptedStageSevenRows'):source.index('TEST(RetainedComputeAdmission, FmaSelectsOnlyNativeWideProducts')]
assert tests.count('TEST(')==6
headers=source[:source.index('#ifdef SIRIUS_HAS_RETAINED_COMPUTE')]
products=Path('bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h').read_text()
row=re.search(r'kTransportRowWords = (\d+)',products).group(1)
subset=headers+'namespace sirius::backend::retained_program { inline constexpr std::size_t kTransportRowWords = '+row+'; }\nnamespace {\nusing namespace sirius::backend;\n'+cls+'\n'+tests+'\n}\n'
(work/'host_protocol_tests.cpp').write_text(subset)
sha=lambda x:hashlib.sha256(x.encode()).hexdigest()
(work/'subset-join.json').write_text(json.dumps({'source_sha256':sha(source),'class_sha256':sha(cls),'tests_sha256':sha(tests),'tests':6,'transport_row_words':int(row),'subset_sha256':sha(subset),'scope':'Exact current source probe/helper and six actual test bodies. External protocol model only; no Vulkan/scientific claim.'},indent=2)+'\n')
