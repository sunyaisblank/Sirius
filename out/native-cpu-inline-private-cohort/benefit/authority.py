"""Pin exercised AO and original AN benefit authorities; bind current AS paths."""
from pathlib import Path
import hashlib

_AR_SOURCE = Path.cwd() / 'out/native-cpu-packed-math-cohort/benefit'
_AR_PINS = {'authority.py': 'b31465bf412fe93a471500f9c472cf0ad89f0e34cd5f878d4b637ed3b19d5514', 'correctness.py': '12b57e38add9319a91796ce3c3308758161fb26cd9b99641dde6e0715521aa83', 'prepare.py': 'cc9969602139fe40b86e41d1476a9e587052c5eda1783f28b67baebde87b80de', 'readers/decimal_check.py': '37705aa878ac58daa25d6f14588bdc06be0dbec8a825e0aa59103844aa7f3928', 'readers/native_evidence.py': '76dbffcce36e1d4a11161085631e06be1aa0e561046ca4b404201b624a22c453', 'readers/numeric_gates.py': '78c899fde5b7846e28c857323fa68af8fc70e8c6a24097e404a70cbb9acbdb2b', 'readers/read_results.py': 'a65b61300702d2d8f254341b5457d8b5246ce10cfb8534c852540fef6e786a26', 'readers/run_abba.py': 'b5b025f0a6ced48e703dd558a7b9c2405ae1b8ddb3030b3a4743bbae1efe2bc4', 'runtime/owned_guard.py': '19f2ba5d29da31c8be60d3644730c0b634288565f212fcafa872354d0004aaee', 'runtime/run_case.py': 'e451025d1518eb558adce27d9369d1f27ebd4aa9bc03af818856cfb98e5d1086', 'runtime/run_owned.py': '22503c90536d4f93fb1f8c22a3a6c39a9b8237f7a1099013346bafd9a24a907f', 'readers/expected-readback-layout.json': 'a40ce35b744ff05c987833de8337ec6f5445daea9553348cdb9209ad78314392'}
for _name, _digest in _AR_PINS.items():
    assert hashlib.sha256((_AR_SOURCE / _name).read_bytes()).hexdigest() == _digest, _name
AP_PINS = {'authority.py': '2e73c790b4fe500788b3de9cdad0b41d279bb85af6f94a913a39ef45f945916f', 'correctness.py': '7ab6d0815407ef6c6b7a95249b2fdc954d9cd50496ab606bc34585cf1b540e13', 'prepare.py': '0ea6c40fcc4861d362830dc8c67c83550c6cdf7185b112be714b31345238d3b0', 'readers/decimal_check.py': '37705aa878ac58daa25d6f14588bdc06be0dbec8a825e0aa59103844aa7f3928', 'readers/native_evidence.py': '76dbffcce36e1d4a11161085631e06be1aa0e561046ca4b404201b624a22c453', 'readers/numeric_gates.py': '78c899fde5b7846e28c857323fa68af8fc70e8c6a24097e404a70cbb9acbdb2b', 'readers/read_results.py': 'd3225b9e35e538a34cc1f8e04631e18ca2e1ace4aee86d02d7e45cddff5186fd', 'readers/run_abba.py': '1e4c93726d87efb029feaa86a900916e5b573a1576cf1e043d26742549a5d770', 'runtime/owned_guard.py': '19f2ba5d29da31c8be60d3644730c0b634288565f212fcafa872354d0004aaee', 'runtime/run_case.py': '1ed4d000c28cd0cf359d321d17e2733543c431f6d3af865ed8fbc13539045613', 'runtime/run_owned.py': '22503c90536d4f93fb1f8c22a3a6c39a9b8237f7a1099013346bafd9a24a907f', 'readers/expected-readback-layout.json': 'a40ce35b744ff05c987833de8337ec6f5445daea9553348cdb9209ad78314392'}
AP_SOURCE = Path.cwd() / 'out/native-cpu-endpoint-header-cohort/benefit'
for name, digest in AP_PINS.items():
    assert hashlib.sha256((AP_SOURCE / name).read_bytes()).hexdigest() == digest, name

AO_PINS = {'prepare.py': 'b413f741f0fc97afae6d7a018727ecd12bd1873832d302ad810166fa3ed5bda6', 'readers/decimal_check.py': 'd36ffe28581c7577e46bcdd77884fcc8475673d8bff37e1a4855c6008d4c4af5', 'readers/native_evidence.py': '69e565e97b5b4ae034d95ac824ddfc02b43e45bdfa1bb583778f820f806ec3ae', 'readers/numeric_gates.py': '18e13a644ef01f83d8818ce37de76c2a6639474266ec61560b4541db60e5379f', 'readers/read_results.py': '1ef811de347f5002dd0a8fc75fc445c458fbb639de3b6e23b94699bebdd27723', 'readers/run_abba.py': '0f77beb5d78cb70eb893f9bef65823dd0b87a7a8405b072e2d2119910594d121', 'runtime/run_case.py': 'd8b695f923ee7eb0cb6c15e3c1c1f4960ed4ed9f911144008d40ef52df113b51', 'runtime/run_owned.py': '95127cdfbdae57cb6b08f75b90a66d45c06f4acded269064f2645e791b619347', 'runtime/owned_guard.py': '5a578e7ed136a69ac15ccb6fc14c2de71d17d2e09ed93fbe2e028ad248ab3042'}
AN_PINS = {'prepare.py': '2384f8db318cc3d5623cf65066ff609a4eb5a0485f87f2704fab60441a44000c', 'readers/decimal_check.py': 'c1861053cb6216ef91f13852b3f3447d03a8d665ba30f3f0c1bf803502bd6ac3', 'readers/native_evidence.py': 'bf50b12e63d42c0715007e173e25ea2abbf01f71030b9d62099be034d59e40e8', 'readers/numeric_gates.py': '5317872eae1fccf0ba9a04f552f3c303e309fda2f0279e573ef1137d8772ca12', 'readers/read_results.py': '3f3f1a7e3703a084460ad077777ebfd4defa5a1290fe9586ef5f664765dfe07b', 'readers/run_abba.py': 'a5b93d12093307cb8f80150c6b74ad2364455f4a697e871417f49ee695ba90a9', 'runtime/run_case.py': '286d5cfa6c9644e374b147595d4ad4d2b31d9c3b68e0350749dbb1fcf7364912', 'runtime/run_owned.py': '95127cdfbdae57cb6b08f75b90a66d45c06f4acded269064f2645e791b619347', 'runtime/owned_guard.py': '5a578e7ed136a69ac15ccb6fc14c2de71d17d2e09ed93fbe2e028ad248ab3042'}


def source(path, digest):
    data = path.read_bytes()
    assert hashlib.sha256(data).hexdigest() == digest, path
    return data.decode()


def load(name):
    text = source(Path.cwd() / 'out/native-cpu-v4-cohort/benefit' / name, AO_PINS[name])
    return text.replace('out/native-cpu-v4-cohort', 'out/native-cpu-inline-private-cohort').replace('out/native-cpu-v4-provider/provider', 'out/native-cpu-inline-private-provider/provider')


def an(name):
    return source(Path.cwd() / 'out/native-cpu-local-identities-cohort/benefit' / name, AN_PINS[name])


def authority_files():
    root = Path.cwd()
    return ([root / 'out/native-cpu-v4-cohort/benefit' / name for name in AO_PINS] +
            [root / 'out/native-cpu-local-identities-cohort/benefit' / name for name in AN_PINS] +
            [root / 'out/native-cpu-active-prefix/benefit' / name for name in AN_PINS] +
            [AP_SOURCE / name for name in AP_PINS])


_previous_authority_files = authority_files
def authority_files():
    return _previous_authority_files() + [_AR_SOURCE / name for name in _AR_PINS]


def prepare_source(text):
    anchor = 'assert all((root / name).is_file() for name in inputs)'
    assert text.count(anchor) == 1
    text = text.replace(anchor,
        'inputs |= {str(path.relative_to(root)) for path in authority_files()}\n'
        'inputs |= {str(path.relative_to(root)) for path in correctness_inputs}\n' + anchor, 1)
    text = text.replace('current H host/AN', 'current H host/AS')
    return text.replace('remaining7/frame only after22PASS.',
        'Original seven controls and the final all8 receipt remain required after fresh22+Decimal; no frame/scientific/release acceptance transfers.')


def reader_scope(text):
    before = 'Remaining7 and cold frame required after22 PASS.'
    assert text.count(before) == 1
    return text.replace(before, 'Current AS preserves AP first7 and adds compiler-private math banks with forced opcode inlining/privateD zeroing and unchanged reviewed zero-residual arithmetic to the same-eight-lane raw-v4 path on the same genuine H host. Original seven controls and final all8 receipt remain required after fresh22+Decimal. Cold-frame/scientific/release gates remain separate.', 1)
