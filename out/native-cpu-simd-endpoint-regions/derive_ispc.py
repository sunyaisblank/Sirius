"""Split only fixed Endpoint induction regions; reuse the sealed K derivation."""
from pathlib import Path
import hashlib
import json
import re
import sys

root = Path.cwd()
work = Path(__file__).resolve().parent
k_generator = root / 'out/native-cpu-simd-structural-uniform/derive_ispc.py'
k_receipt = k_generator.parent / 'ispc-derivation.json'


def sealed(path, sha256):
    data = path.read_bytes()
    assert hashlib.sha256(data).hexdigest() == sha256, path
    return data


k_source = sealed(k_generator, '50495ab7fb902f4750883d6afc31f339c8a945400851894853fbd05a52978e59').decode()
k_receipt_bytes = sealed(k_receipt, '8f7c959c680ef1d982acf1837c148e1ed8332aa13c845a4c38f6fda76c3c0a72')
k_derivation = json.loads(k_receipt_bytes)


def endpoint_regions(text, slug):
    if slug != 'endpoint':
        return text, {'changed': False, 'reused_exact_K_source': True}
    before = text

    def function(name):
        match, = re.finditer(r'^static[^\n]*\b' + name + r'\([^\n]+\)\n\{', before, re.M)
        end = end_of_block(before, match.end() - 1)
        return match.start(), end, before[match.start():end]

    def strip_lines(body):
        return re.sub(r'^\s*#line[^\n]*\n', '', body, flags=re.M)

    def qualify(body, names, allowed):
        for name in names:
            for rhs in re.findall(r'\b' + name + r'\s*=(?!=)\s*([^;\n]+);', body):
                expr = re.sub(r'\b(?:0x[0-9a-fA-F]+|\d+)(?:U|UL|ULL)?\b', '', rhs)
                identifiers = set(re.findall(r'\b[A-Za-z_]\w*\b', expr))
                assert identifiers <= set(allowed), (name, rhs, identifiers - set(allowed))
            body, count = re.subn(r'\buint32_t ' + name + r'\b', 'uniform uint32_t ' + name, body)
            assert count == 1, (name, count)
        return body

    def unqualify(body, names):
        for name in names:
            body, count = re.subn(r'\buniform uint32_t ' + name + r'\b', 'uint32_t ' + name, body)
            assert count == 1
        return body

    mappings = {'ReadStored_0': 'ReadStoredEndpointRegion_0',
                'WriteStored_0': 'WriteStoredEndpointRegion_0',
                'Tangent_0': 'TangentEndpointRegion_0'}
    helper_parameters = {'ReadStored_0': ('base_1',), 'WriteStored_0': ('base_0',),
                         'Tangent_0': ('i_1',)}
    additions, edits = [], []
    for old, selected in helper_parameters.items():
        _, end, original = function(old)
        clone = qualify(original, selected, selected)
        clone = clone.replace(old + '(', mappings[old] + '(', 1)
        if old == 'Tangent_0':
            clone, count = re.subn(r'\bReadStored_0(?=\()', mappings['ReadStored_0'], clone)
            assert count == 1
        restored = clone.replace(mappings[old] + '(', old + '(', 1)
        if old == 'Tangent_0':
            restored = restored.replace(mappings['ReadStored_0'] + '(', 'ReadStored_0(', 1)
        assert unqualify(restored, selected) == original
        addition = '\n\n// Fixed Endpoint induction columns; original generic helper retained.\n' + clone
        additions.append(addition)
        edits.append((end, end, addition))

    start, end, owned = function('ComputeOwnedRow_0')
    stripped = strip_lines(owned)
    original_bounds = re.findall(r'if\((\w+) < (\w+)\)\s*\{\s*\}\s*else\s*\{\s*break;\s*\}', stripped)
    assert original_bounds == [('i_11', '45U'), ('i_11', '53U'), ('phase_index_0', '3U'),
        ('i_11', '40U'), ('axis_0', '4U'), ('i_11', '4U'), ('i_11', '4U'), ('axis_0', '4U')]
    assert len(re.findall(r'\bbreak;', stripped)) == len(original_bounds)
    assert not re.search(r'\bcontinue\s*;', stripped)
    assert not re.search(r'&\s*(?:i_11|axis_0)\b', owned)
    revised, regions = owned, []

    def split(counter, initial, bound, fresh):
        nonlocal revised
        matches = []
        for match in re.finditer(r'(?m)^(?P<indent> +)for\(;;\)\n(?P=indent)\{', revised):
            finish = end_of_block(revised, match.end() - 1)
            loop = revised[match.start():finish]
            normalized = strip_lines(loop)
            first = re.match(r'\s*for\(;;\)\s*\{\s*if\((\w+) < (\w+)\)', normalized)
            reset = re.search(r'(?m)^ *' + counter + r' = ' + initial +
                              r';(?:(?:\s+)|(?:#line[^\n]*\n))*\Z', revised[:match.start()])
            if first and first.groups() == (counter, bound) and reset:
                matches.append((match, finish, loop))
        # The two distinct i<4 regions are processed in original source order.
        assert len(matches) == (2 if fresh == 'endpoint_projection_column' else 1), (fresh, len(matches))
        match, finish, original = matches[0]
        prefix = revised[:match.start()]
        reset = re.search(r'(?m)^ *' + counter + r' = ' + initial + r';(?:(?:\s+)|(?:#line[^\n]*\n))*\Z', prefix)
        assert reset, (counter, initial, bound)
        writes = re.findall(r'\b' + counter + r'\s*=(?!=)\s*([^;\n]+);', original)
        expected = ['axis_1'] if fresh == 'endpoint_flat_axis' else [counter + ' + 1U']
        assert writes == expected, (fresh, writes)
        if fresh == 'endpoint_flat_axis':
            assert original.count('uint32_t axis_1 = axis_0 + 1U;') == 1
        # No helper receives the old variable by address; only this loop and
        # its nested fixed loops can observe the local induction value.
        renamed = re.sub(r'\b' + counter + r'\b', fresh, original)
        indent = match['indent']
        replacement = (indent + '{\n' + indent + '    uniform uint32_t ' + fresh + ' = ' + initial + ';\n'
                       + renamed + '\n' + indent + '    ' + counter + ' = ' + fresh + ';\n' + indent + '}')
        revised = revised[:match.start()] + replacement + revised[finish:]
        regions.append({'original': original, 'replacement': replacement, 'counter': counter,
                        'initial': initial, 'bound': bound, 'fresh': fresh,
                        'original_reset_retained': True, 'masked_final_value_copy_back': True})

    # Inner loops first: outer-loop inverse then recovers these exact regions.
    split('axis_0', '1U', '4U', 'endpoint_flat_axis')
    split('axis_0', '0U', '4U', 'endpoint_matrix_axis')
    split('i_11', '0U', '45U', 'endpoint_input_column')
    split('i_11', '45U', '53U', 'endpoint_auxiliary_column')
    split('i_11', '0U', '40U', 'endpoint_flat_column')
    split('i_11', '0U', '4U', 'endpoint_projection_column')
    split('i_11', '0U', '4U', 'endpoint_matrix_column')

    aliases = ('_S331', '_S332', '_S334', '_S349', '_S365', 'axis_1', '_S385', '_S387',
               '_S389', '_S392', '_S395', '_S398')
    uniform_columns = set(aliases) | {r['fresh'] for r in regions} | {'rowBase_11'}
    revised = qualify(revised, aliases, uniform_columns)
    call_counts = {}
    for old, new in mappings.items():
        count = 0
        for match in re.finditer(r'\b' + old + r'\(', revised):
            opening = match.end() - 1
            closing = close_delimiter(revised, opening, '(', ')')
            arguments = revised[opening + 1:closing]
            column = arguments.split(',')[1 if old == 'Tangent_0' else 0]
            expr = re.sub(r'\b\d+U\b', '', column)
            identifiers = set(re.findall(r'\b[A-Za-z_]\w*\b', expr))
            assert identifiers <= uniform_columns, (old, column, identifiers - uniform_columns)
            count += 1
        assert count == {'ReadStored_0': 9, 'WriteStored_0': 8, 'Tangent_0': 2}[old]
        revised, changed = re.subn(r'\b' + old + r'(?=\()', new, revised)
        assert changed == count
        call_counts[old] = count
    restored = revised
    for old, new in mappings.items():
        restored = re.sub(r'\b' + new + r'(?=\()', old, restored)
    restored = unqualify(restored, aliases)
    for region in reversed(regions):
        assert restored.count(region['replacement']) == 1
        restored = restored.replace(region['replacement'], region['original'], 1)
    assert restored == owned, 'Endpoint owned-body inverse must restore every K byte'
    edits.append((start, end, revised))
    for start, end, replacement in sorted(edits, reverse=True):
        text = text[:start] + replacement + text[end:]
    restored = text.replace(revised, owned, 1)
    for addition in additions:
        assert restored.count(addition) == 1
        restored = restored.replace(addition, '', 1)
    assert restored == before, 'Endpoint full-source inverse must restore every K byte'
    return text, {'changed': True, 'regions': [
        {key: value for key, value in region.items() if key not in ('original', 'replacement')}
        for region in regions], 'uniform_output_aliases': list(aliases),
        'owned_helper_replacements': call_counts,
        'helper_clones': [{'function': mappings[name], 'uniform_parameters': list(selected)}
                          for name, selected in helper_parameters.items()],
        'generic_helpers_retained': True, 'full_inverse_restores_K_source_bytes': True,
        'all_loop_breaks_are_original_fixed_bounds': True,
        'projection_success_component_and_handoff_remain_varying': True,
        'phase_counter_and_evaluator_unchanged': True,
        'early_return_domain': 'A numerical return retires that row from the owned call. Every continuing row in a lexical region has the same constant initial value, fixed bound and increment. Final counter copy-back occurs under the same surviving mask; old counters are not aliased. Per-row physical values and predicates remain varying.'}


def replace_once(old, new):
    global k_source
    assert k_source.count(old) == 1, old[:100]
    k_source = k_source.replace(old, new, 1)


replace_once("destination = root / 'bin/linux-gcc/native-cpu-simd-structural-uniform/ispc-source'",
             "destination = root / 'bin/linux-gcc/native-cpu-simd-endpoint-regions/ispc-source'")
replace_once('    text, structural = structural_uniform(text, slug)\n',
             '    text, structural = structural_uniform(text, slug)\n'
             '    k_text = text\n    text, region_proof = endpoint_regions(text, slug)\n')
replace_once("    derived = PRELUDE + text + entry\n    if slug == 'endpoint':\n        assert pin(derived.encode()) == j_derivation['stages'][slug]['derived']\n",
             "    assert pin((PRELUDE + k_text + entry).encode()) == k_derivation['stages'][slug]['derived']\n"
             "    derived = PRELUDE + text + entry\n    if slug != 'endpoint':\n"
             "        assert pin(derived.encode()) == k_derivation['stages'][slug]['derived']\n")
replace_once("        'structural_uniformity': structural,\n",
             "        'structural_uniformity': structural,\n        'endpoint_regions': region_proof,\n")
replace_once("{'instructions': total, 'stages': records,\n",
             "{'instructions': total, 'stages': records,\n"
             "    'derivation_basis': {'generator_path': str(k_generator.relative_to(root)), 'generator': pin(k_generator.read_bytes()),\n"
             "                         'receipt_path': str(k_receipt.relative_to(root)), 'receipt': pin(k_receipt_bytes)},\n")
exec(compile(k_source, str(k_generator), 'exec'), globals())
