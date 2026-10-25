"""Preserve all original C8 numerical/module checks; add actual field route proof."""
from pathlib import Path
basis=Path.cwd()/'out/native-cpu-active-prefix/read_result.py'
source=basis.read_text().replace('_try_packet','_rows_try_packet').replace('providerB','providerT')
before="endpoint, = [row for row in observation['calls']"
assert source.count(before)==1
source=source.replace(before,"""namespace={}
proof_source=work/'route_proof.py'
exec(compile(proof_source.read_text(),str(proof_source),'exec'),namespace)
field_route_proof=namespace['check'](root,action,observation)
"""+before)
source=source.replace("'native_first_packets_handled': 3,","'native_first_packets_handled': 3, 'field_route_proof': field_route_proof,")
source=source.replace('patch changes objectpaths/cacheUUID only','patch replaces row SIMD by field SIMD and original logical row tasks; unchanged feature/property implementation')
exec(compile(source,str(basis),'exec'))
