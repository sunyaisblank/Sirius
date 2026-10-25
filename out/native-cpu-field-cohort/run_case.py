"""Use the original one-C8 runner; require each new observed route's return."""
from pathlib import Path
basis=Path.cwd()/'out/native-cpu-active-prefix/run_case.py'
source=basis.read_text().replace('_try_packet','_rows_try_packet')
before="assert (folder / 'gtest.xml').is_file()"
assert source.count(before)==1
source=source.replace(before,"""for category in ('preflight_calls','field_calls','layer_calls'):
    assert len(record[category])==3
    assert all(row['returned'] and row['handled'] for row in record[category])
capability,=record['capability_calls']
assert capability['returned'] and capability['handled']
"""+before)
exec(compile(source,str(basis),'exec'))
