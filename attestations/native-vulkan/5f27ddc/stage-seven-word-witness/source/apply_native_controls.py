"""Apply the original native controls to one compiler assembly; no device."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import ast
import re

work = Path(__file__).resolve().parent
assert len(sys.argv) == 2 and sys.argv[1] in ['transport', 'endpoint']
stage = sys.argv[1]
path = Path('scripts/build-retained-kernels.py')
tree = ast.parse(path.read_text(), str(path))
functions = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == 'native_controls']
assert len(functions) == 1
namespace = {'re': re}
exec(compile(ast.Module(body=functions, type_ignores=[]), str(path), 'exec'), namespace)
source = work / 'kernels' / ('retained_' + stage + '.compiler.spvasm')
destination = work / 'kernels' / ('retained_' + stage + '.spvasm')
assert not destination.exists()
text = namespace['native_controls'](source.read_text(), False)
assert 'OpReadClockKHR' not in text and 'OpCapability ShaderClockKHR' not in text
assert 'DescriptorSet 0' in text and 'Binding 2' in text
destination.write_text(text)
