"""Original complete retained import/include closure; historical/current scope stays distinct."""
import ast,re
from pathlib import Path

def kernel_source_profile(P,R,W,record,content):
    old = Path(P['original_stock_work']) / 'stock/src/sirius/kernels'
    private = W / 'stock/src/sirius/kernels'
    canonical = R / 'src/sirius/kernels'
    names = sorted(p.name for p in old.iterdir())
    assert len(names) == 33 and names == sorted(p.name for p in private.iterdir())
    assert names == sorted(p.name for p in (W / 'candidate/src/sirius/kernels').iterdir())
    rows, deltas = {}, []
    for name in names:
        stock = private / name; current = canonical / name
        assert stock.read_bytes() == (old / name).read_bytes(), name
        if name != 'portable_binary32.h':
            assert stock.read_bytes() == (W / 'candidate/src/sirius/kernels' / name).read_bytes(), name
        different = stock.read_bytes() != current.read_bytes()
        if different:
            deltas.append(name)
            expected = P['known_outside_retained_source_deltas'][name]
            assert content(stock) == expected['original_private'] and content(current) == expected['canonical4498'], name
        rows[name] = {'original':record(old/name), 'private':record(stock), 'canonical4498':record(current), 'current_equal':not different}
    assert sorted(deltas) == sorted(P['known_outside_retained_source_deltas']) == ['gr_ellis_topology.slang','gr_integrator.slang','parity_probe.slang','trace.slang']
    # Literal Slang import/include closure only. Every directive is recognized;
    # macro/dynamic includes, include_next and alternate search names refuse.
    seen, pending, generated = set(), list(P['retained_source_roots']), set()
    while pending:
        name = pending.pop()
        if name in seen: continue
        seen.add(name); assert name in names
        text = (private/name).read_text()
        if name.endswith('.py'):
            tree = ast.parse(text)
            imports = []
            for node in ast.walk(tree):
                if isinstance(node,ast.Import): imports.extend(x.name for x in node.names)
                elif isinstance(node,ast.ImportFrom):
                    assert node.level == 0; imports.append(node.module)
                elif isinstance(node,ast.Call):
                    assert not (isinstance(node.func,ast.Name) and node.func.id in ('open','eval','exec','compile','__import__','getattr'))
                    assert not (isinstance(node.func,ast.Attribute) and node.func.attr in ('read_text','read_bytes','open','read','write_text','write_bytes','exec_module','import_module'))
            assert sorted(imports) == ['fractions','struct']
            continue
        text = re.sub(r'/\*.*?\*/|//[^\n]*','',text,flags=re.S)
        for line in text.splitlines():
            if re.search(r'\bimport\b',line):
                match = re.fullmatch(r'\s*import (\w+);\s*',line); assert match, line
                dependency = match[1]+'.slang'
                if dependency == 'retained_owner_layout.slang':
                    assert name == 'retained_owner.slang'; generated.add(dependency)
                else: pending.append(dependency)
            if re.search(r'^\s*#\s*(?:include|include_next|import)\b',line):
                if name == 'portable_binary32.h' and line in ('#include <bit>','#include <cstdint>'):
                    cpp_branch = text.split('#if defined(__cplusplus)',1)[1].split('#else',1)[0]
                    assert line in cpp_branch.splitlines()
                    continue  # Fixed C++-only headers, inactive for retained Slang.
                match = re.fullmatch(r'\s*#include "([A-Za-z0-9_.]+)"\s*',line); assert match, line
                pending.append(match[1])
    assert sorted(seen) == P['retained_source_closure'] and len(seen) == 14
    assert generated == {'retained_owner_layout.slang'} and not set(deltas)&seen
    assert all(rows[name]['current_equal'] for name in seen)
    builder = (W / 'stock/scripts/build-retained-kernels.py').read_text()
    assert builder == (R / 'scripts/build-retained-kernels.py').read_text()
    assert builder == (Path(P['original_stock_work']) / 'stock/scripts/build-retained-kernels.py').read_text()
    assert builder.count('spec_from_file_location("retained_program", source / "retained_program.py")') == 1
    assert 'owner_include = ["-I", str(destination.parent)] if stage_seven_owner else []' in builder
    assert '[compiler, str(source), *definitions, "-I", str(source.parent), *owner_include, optimization,' in builder
    # Both search roots are fixed by this source: private source directory and
    # owned generated owner-layout directory. No CLI source/include argument.
    cmake = (R / 'src/sirius/backend/CMakeLists.txt').read_text()
    deps = cmake.split('DEPENDS "${CMAKE_SOURCE_DIR}/scripts/build-retained-kernels.py"',1)[1].split('COMMENT "Building bounded retained camera and joint transport stages"',1)[0]
    assert sorted(re.findall(r'"\$\{SIRIUS_KERNEL_DIR\}/([^"/]+)"',deps)) == sorted(seen)
    layout = lambda source: source.split('def generate_owner_layout(',1)[1].split('\ndef compile_shader(',1)[0]
    assert layout(builder) == layout((W / 'candidate/scripts/build-retained-kernels.py').read_text())
    return {'historical33_original_private_sources':rows, 'complete_retained_import_include_closure':sorted(seen),
            'generated_owner_layout':sorted(generated),'inactive_cpp_only_standard_headers':['bit','cstdint'],'all14_original_stock_consumed_source_files_equal_canonical4498':True,
            'current_equal_source_count':len(names)-len(deltas), 'exact4_current_deltas_outside_retained_closure':deltas,
            'builder':record(W/'stock/scripts/build-retained-kernels.py'),'CMake_dependency_authority':record(R/'src/sirius/backend/CMakeLists.txt'),
            'scope':'Historical4bcd complete stock sourcebase preserved, complete retained original-stock actual-read14 closure equal current4498; candidate header is the separately pinned accepted private insertion. Four unrelated floor-repaired trace/probe files are frozen current deltas outside that closure. No current full trace/science qualification.'}
