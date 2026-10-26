"""Pin exercised AO and original AN benefit authorities; bind current AP paths."""
from pathlib import Path
import hashlib
AO_PINS = {'prepare.py': 'b413f741f0fc97afae6d7a018727ecd12bd1873832d302ad810166fa3ed5bda6', 'readers/decimal_check.py': 'd36ffe28581c7577e46bcdd77884fcc8475673d8bff37e1a4855c6008d4c4af5', 'readers/native_evidence.py': '69e565e97b5b4ae034d95ac824ddfc02b43e45bdfa1bb583778f820f806ec3ae', 'readers/numeric_gates.py': '18e13a644ef01f83d8818ce37de76c2a6639474266ec61560b4541db60e5379f', 'readers/read_results.py': '1ef811de347f5002dd0a8fc75fc445c458fbb639de3b6e23b94699bebdd27723', 'readers/run_abba.py': '0f77beb5d78cb70eb893f9bef65823dd0b87a7a8405b072e2d2119910594d121', 'runtime/run_case.py': 'd8b695f923ee7eb0cb6c15e3c1c1f4960ed4ed9f911144008d40ef52df113b51', 'runtime/run_owned.py': '95127cdfbdae57cb6b08f75b90a66d45c06f4acded269064f2645e791b619347', 'runtime/owned_guard.py': '5a578e7ed136a69ac15ccb6fc14c2de71d17d2e09ed93fbe2e028ad248ab3042'}
AN_PINS = {'prepare.py': '2384f8db318cc3d5623cf65066ff609a4eb5a0485f87f2704fab60441a44000c', 'readers/decimal_check.py': 'c1861053cb6216ef91f13852b3f3447d03a8d665ba30f3f0c1bf803502bd6ac3', 'readers/native_evidence.py': 'bf50b12e63d42c0715007e173e25ea2abbf01f71030b9d62099be034d59e40e8', 'readers/numeric_gates.py': '5317872eae1fccf0ba9a04f552f3c303e309fda2f0279e573ef1137d8772ca12', 'readers/read_results.py': '3f3f1a7e3703a084460ad077777ebfd4defa5a1290fe9586ef5f664765dfe07b', 'readers/run_abba.py': 'a5b93d12093307cb8f80150c6b74ad2364455f4a697e871417f49ee695ba90a9', 'runtime/run_case.py': '286d5cfa6c9644e374b147595d4ad4d2b31d9c3b68e0350749dbb1fcf7364912', 'runtime/run_owned.py': '95127cdfbdae57cb6b08f75b90a66d45c06f4acded269064f2645e791b619347', 'runtime/owned_guard.py': '5a578e7ed136a69ac15ccb6fc14c2de71d17d2e09ed93fbe2e028ad248ab3042'}


def source(path, digest):
    data = path.read_bytes()
    assert hashlib.sha256(data).hexdigest() == digest, path
    return data.decode()


def load(name):
    text = source(Path.cwd() / 'out/native-cpu-v4-cohort/benefit' / name, AO_PINS[name])
    return text.replace('out/native-cpu-v4-cohort', 'out/native-cpu-endpoint-header-cohort').replace('out/native-cpu-v4-provider/provider', 'out/native-cpu-endpoint-header-provider/provider')


def an(name):
    return source(Path.cwd() / 'out/native-cpu-local-identities-cohort/benefit' / name, AN_PINS[name])


def authority_files():
    root = Path.cwd()
    return ([root / 'out/native-cpu-v4-cohort/benefit' / name for name in AO_PINS] +
            [root / 'out/native-cpu-local-identities-cohort/benefit' / name for name in AN_PINS] +
            [root / 'out/native-cpu-active-prefix/benefit' / name for name in AN_PINS])


def prepare_source(text):
    anchor = 'assert all((root / name).is_file() for name in inputs)'
    assert text.count(anchor) == 1
    text = text.replace(anchor,
        'inputs |= {str(path.relative_to(root)) for path in authority_files()}\n'
        'inputs |= {str(path.relative_to(root)) for path in correctness_inputs}\n' + anchor, 1)
    text = text.replace('current H host/AN', 'current H host/AP')
    return text.replace('remaining7/frame only after22PASS.',
        'Original seven controls completed before this trial; final all8 acceptance still waits for fresh22+Decimal and no frame/scientific/release acceptance transfers.')


def reader_scope(text):
    before = 'Remaining7 and cold frame required after22 PASS.'
    assert text.count(before) == 1
    return text.replace(before, 'Current AP Endpoint admission composition uses six unchanged AO objects and the same genuine H host. Original seven controls completed before this trial; final all8 receipt still requires fresh22+Decimal. A cold frame and scientific/release gates remain separate.', 1)
