"""Independent Decimal algebra on exact bytes from one accepted failed Kerr capture."""
from decimal import Decimal, localcontext
import hashlib
import json
import struct
import sys


if not __debug__:
    raise RuntimeError('Assertions must be enabled')

SOURCE = '02ee51f85e45e8342923bae1d9f26d50f1d39e6d'


def raw_hex(value, size):
    assert type(value) is str and len(value) == 2 * size
    result = bytes.fromhex(value)
    assert len(result) == size and result.hex() == value
    return result


def exact(value):
    result = Decimal.from_float(value)
    assert result.is_finite()
    return result


def record(raw):
    values = struct.unpack_from('<60f', raw)
    vectors = [values[index:index + 4] for index in range(0, 52, 4)]
    names = ['x', 'u', 'X0', 'X1', 'X2', 'X3', 'V0', 'V1', 'V2', 'V3',
             'compX', 'compU', 'integration']
    return {
        'vectors': {name: [str(exact(value)) for value in vector]
                    for name, vector in zip(names, vectors)},
        'source': [str(exact(value)) for value in values[52:56]],
        'volume': [str(exact(value)) for value in values[56:60]],
        'control': list(struct.unpack_from('<4I', raw, 240)),
        'identity': list(struct.unpack_from('<4I', raw, 256)),
    }


def inverse_and_determinant(matrix):
    # Generic pivoted elimination; no Kerr inverse or null-cone assumption.
    augmented = [list(row) + [Decimal(int(i == j)) for j in range(4)]
                 for i, row in enumerate(matrix)]
    determinant = Decimal(1)
    for column in range(4):
        pivot_row = max(range(column, 4), key=lambda row: abs(augmented[row][column]))
        assert augmented[pivot_row][column] != 0, 'Singular independent metric'
        if pivot_row != column:
            augmented[column], augmented[pivot_row] = augmented[pivot_row], augmented[column]
            determinant = -determinant
        pivot = augmented[column][column]
        determinant *= pivot
        augmented[column] = [value / pivot for value in augmented[column]]
        for row in range(4):
            if row != column:
                factor = augmented[row][column]
                augmented[row] = [value - factor * other
                                  for value, other in zip(augmented[row], augmented[column])]
    inverse = [row[4:] for row in augmented]
    residual = max(abs(sum(matrix[i][k] * inverse[k][j] for k in range(4)) - int(i == j))
                   for i in range(4) for j in range(4))
    return determinant, residual


def evaluate(raw, packet, precision):
    # Independent defining Cartesian Kerr-Schild geometry, as in the separated
    # reference's At(). No product imports, Christoffels, projector or controller.
    with localcontext() as context:
        context.prec = precision
        values = [exact(value) for value in struct.unpack_from('<60f', raw)]
        mass, spin, charge = (exact(packet[index]) for index in (3, 4, 5))
        x, y, z = values[1:4]
        tangent = values[4:8]
        rho2 = x*x + y*y + z*z
        reduced = rho2 - spin*spin
        discriminant = reduced*reduced + 4*spin*spin*z*z
        root = discriminant.sqrt()
        # Positive root of r^4-(rho^2-a^2)r^2-a^2*z^2=0; rationalize
        # the negative-reduced branch instead of losing its positive root.
        radius2 = ((reduced + root) / 2 if reduced >= 0
                   else 2*spin*spin*z*z / (root - reduced))
        result = {'decimal_precision': precision,
                  'radius_squared': str(radius2),
                  'quartic_residual': str(radius2*radius2 - reduced*radius2 - spin*spin*z*z),
                  'stored_adaptive_step': str(values[48]),
                  'parameter_minimum_step': str(exact(packet[23])),
                  'parameter_maximum_step': str(exact(packet[24])),
                  'parameter_step_scale': str(exact(packet[22])),
                  'stored_adaptive_step_equals_floor': values[48] == exact(packet[23]),
                  'geometry_evaluable': radius2 > 0}
        if radius2 <= 0:
            result['geometry_issue'] = 'No positive Kerr radial root at this recorded event'
            return result
        radius = radius2.sqrt()
        sigma = radius2 + spin*spin*z*z/radius2
        coefficient = (2*mass*radius - charge*charge) / sigma
        ell = [Decimal(1), (radius*x + spin*y)/(radius2 + spin*spin),
               (radius*y - spin*x)/(radius2 + spin*spin), z/radius]
        eta = [Decimal(-1), Decimal(1), Decimal(1), Decimal(1)]
        metric = [[(eta[i] if i == j else Decimal(0)) + coefficient*ell[i]*ell[j]
                   for j in range(4)] for i in range(4)]
        terms = [metric[i][j]*tangent[i]*tangent[j] for i in range(4) for j in range(4)]
        contraction = sum(terms)
        absolute_terms = sum(abs(value) for value in terms)
        momentum = [sum(metric[i][j]*tangent[j] for j in range(4)) for i in range(4)]
        determinant, inverse_residual = inverse_and_determinant(metric)
        horizon_discriminant = mass*mass - spin*spin - charge*charge
        horizon = mass + horizon_discriminant.sqrt() if horizon_discriminant >= 0 else None
        capture = horizon*exact(packet[26]) if horizon is not None else None
        frequency = -tangent[0]/(1 + coefficient).sqrt() if 1 + coefficient > 0 else None
        result.update({
            'radius': str(radius), 'sigma': str(sigma), 'metric_h_coefficient': str(coefficient),
            'ell_covector': [str(value) for value in ell],
            'ell_minkowski_norm': str(sum(eta[i]*ell[i]*ell[i] for i in range(4))),
            'metric_covariant': [[str(value) for value in row] for row in metric],
            'metric_determinant': str(determinant),
            'metric_determinant_plus_one': str(determinant + 1),
            'metric_inverse_max_identity_residual': str(inverse_residual),
            'gkk': str(contraction), 'absolute_metric_tangent_term_sum': str(absolute_terms),
            'normalized_null_residual': str(abs(contraction)/absolute_terms) if absolute_terms > 0 else None,
            'covariant_momentum': [str(value) for value in momentum],
            'energy_minus_p_t': str(-momentum[0]),
            'axial_angular_momentum': str(x*momentum[2] - y*momentum[1]),
            'independent_outer_horizon_radius': str(horizon) if horizon is not None else None,
            'capture_factor': str(exact(packet[26])),
            'independent_buffered_capture_radius': str(capture) if capture is not None else None,
            'radius_minus_independent_capture_radius': str(radius - capture) if capture is not None else None,
            'reference_past_frequency': str(frequency) if frequency is not None else None,
            'reference_past_frequency_positive': frequency > 0 if frequency is not None else None,
        })
        return result


def main():
    assert len(sys.argv) == 3, 'usage: analyse_state.py exact-bridge.json new-output.json'
    input_path, output_path = sys.argv[1:]
    assert input_path != output_path
    with open(input_path, 'rb') as stream:
        original = stream.read()
    bridge = json.loads(original)
    assert bridge['source_revision'] == bridge['runner_source_revision'] == SOURCE
    assert bridge['kind'] == 'kerr-image' and bridge['accepted_owner_execution'] is True
    assert bridge['completed'] is True and bridge['expected_completed_failure_captured'] is True
    assert bridge['numerical_test_passed'] is False
    assert type(bridge['returncode']) is int and bridge['returncode'] == 1
    assert type(bridge['expected_child_exit']) is int and bridge['expected_child_exit'] == 1
    assert bridge['whole_inputs_unchanged'] is True and bridge['linux_child_birth_absent'] is True
    assert not bridge['cleanup_errors'] and not bridge['signals']
    assert not bridge.get('error') and not bridge.get('postcondition_error')
    capture = bridge['failure_context']
    assert capture['schema'] == 'sirius-native-windows-kerr-failure-context-v1'
    assert capture['numerical_test_passed'] is False and capture['radiance_publication_admitted'] is False
    raw_context = capture['raw_context_json']
    failure = json.loads(raw_context)
    assert failure['schema'] == 'sirius-trace-continuation-failure-v1'
    assert failure['scalar_bytes'] == failure['packet_scalar_bytes'] == 4
    assert failure['record_bytes'] == 272 and failure['byte_order'] == 'little'
    assert failure['action'] == 1 and failure['active'] == 64
    packet_raw = raw_hex(failure['packet_hex'], 280)
    previous = raw_hex(failure['previous_record_hex'], 272)
    current = raw_hex(failure['current_record_hex'], 272)
    packet = struct.unpack('<70f', packet_raw)
    packet_words = struct.unpack('<70I', packet_raw)
    assert packet[2] == packet[5] == packet[6] == packet[46] == 0.0
    assert packet[3] == 1.0 and 0 < packet[4] < packet[3]
    assert packet[0] == packet[1] == 64.0 and packet[21] == 2000.0 and packet[68] == 1.0
    assert packet[23] > 0 and packet[24] >= packet[23]
    before, after = record(previous), record(current)
    identity = failure['expected_identity']
    assert before['identity'] == after['identity'] == identity
    assert identity[2] == 0 and identity[3] == packet_words[69] > 0
    assert after['control'][0] == 2 and after['control'][1] == 5
    rectangle = failure['active_rectangle']
    assert list(packet[31:35]) == rectangle
    assert rectangle[0] == 0 and 0 <= rectangle[1] < 64 and rectangle[2:] == [64, 1]
    assert type(failure['row']) is int and 0 <= failure['row'] < 64
    assert identity[:2] == [failure['row'], rectangle[1]]
    name_raw = bytes.fromhex(failure['device_name_bytes_hex'])
    assert name_raw.hex() == failure['device_name_bytes_hex']
    assert name_raw.decode('utf-8') == 'AMD Radeon 780M Graphics'
    rollback = {
        'physical_vectors_0_through_11_equal': previous[:192] == current[:192],
        'accepted_affine_length_last_interval_beam_seed_equal': previous[196:208] == current[196:208],
        'source_and_volume_equal': previous[208:240] == current[208:240],
        'accepted_attempt_count_equal': before['control'][3] == after['control'][3],
        'one_attempt_added': after['control'][2] == before['control'][2] + 1,
        'previous_running': before['control'][0] == 1 and before['control'][1] == 0,
    }
    assert rollback == capture['rollback_comparisons']
    calculations = {str(precision): {label: evaluate(raw, packet, precision)
                                    for label, raw in [('previous', previous), ('current', current)]}
                    for precision in (80, 60)}
    differences = {}
    with localcontext() as context:
        context.prec = 80
        for label in ['previous', 'current']:
            differences[label] = {}
            for key in ['radius', 'gkk', 'normalized_null_residual', 'energy_minus_p_t',
                        'axial_angular_momentum', 'metric_determinant']:
                high, low = calculations['80'][label].get(key), calculations['60'][label].get(key)
                differences[label][key] = str(abs(Decimal(high) - Decimal(low))) if high is not None and low is not None else None
    with open(__file__, 'rb') as stream:
        script_sha256 = hashlib.sha256(stream.read()).hexdigest()
    result = {
        'schema': 'sirius-kerr-failure-independent-decimal-state-analysis-v1',
        'analysis_completed': True, 'source_revision': SOURCE, 'numerical_test_passed': False,
        'input': {'path': input_path, 'bytes': len(original), 'sha256': hashlib.sha256(original).hexdigest()},
        'script_sha256': script_sha256, 'raw_context_json': raw_context,
        'packet_exact_decimal': [str(exact(value)) for value in packet],
        'packet_words': list(packet_words), 'previous': before, 'current': after,
        'rollback_comparisons': rollback, 'rollback_is_observation_not_acceptance': True,
        'calculations': calculations, 'absolute_80_minus_60_digit_differences': differences,
        'precision_differences_are_not_certified_error_bounds': True,
        'scope': 'Independent defining Kerr-Schild state algebra only, evaluated separately at80/60 digits from exact binary32 values. No product metric/connection/projector/controller imports, stepping, shader execution or full-ray inference. No conclusion identifies the actual rejected GPU conjunct; trial candidates were not captured. Buffered capture radius is an independent high-precision value, not a claim of bit-identical shader evaluation.',
    }
    with open(input_path, 'rb') as stream:
        assert stream.read() == original, 'Input changed during analysis'
    with open(output_path, 'x', encoding='utf-8') as stream:
        stream.write(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
