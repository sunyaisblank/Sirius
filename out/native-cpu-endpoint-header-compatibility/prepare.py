from authority import load
source = load('prepare.py')
anchor = 'exec(compile(source,str(basis),"exec"))'
assert source.count(anchor) == 1
addition = '''
from authority import authority_files
anchor = "tool_receipt = R / 'out/native-cpu-host-representation/probe/preparation.json'"
assert source.count(anchor) == 1
source = source.replace(anchor, 'for path in authority_files():\\n    inputs[str(path.relative_to(R))] = pin(path)\\n' + anchor, 1)
exec(compile(source, str(basis), "exec"))
'''
source = source.replace(anchor, addition, 1)
exec(compile(source, '/home/astra/.project/Sirius/out/native-cpu-v4-compatibility/prepare.py', 'exec'))
