"""Exact original generated uint32 array parser only."""
import re,struct

def arrays(path):
    result = {}
    for count, name, body in re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};', path.read_text(), re.S):
        words = [int(v) for v in re.findall(r'(\d+)u', body)]
        assert len(words) == int(count) and name not in result
        result[name] = struct.pack('<' + str(len(words)) + 'I', *words)
    return result
