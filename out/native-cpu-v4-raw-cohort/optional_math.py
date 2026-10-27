from rebind import ao
exec(compile(ao('optional_math.py'), '/home/astra/.project/Sirius/out/native-cpu-v4-cohort/optional_math.py', "exec"))

_original_check = check


def check(root, action, observation):
    result = _original_check(root, action, observation)
    assert {'endpoint', 'dopri_phase'} <= set(result['selected_simd_stages'])
    assert result['all_selected_optional_math_observed']
    assert set(result['observed_optional_math_stages']) == set(result['selected_simd_stages'])
    return result
