from pathlib import Path
import json,re

W=Path(__file__).resolve().parent
source=(W/'stock-portable_binary32.h').read_text()
source=source.replace('#include <bit>','#include <bit>\n#include <cmath>')
def function(text,name):
    match=re.search(r'SIRIUS_PB32_INLINE [^\n]+ '+name+r'\(',text)
    assert match,name
    start=match.start(); end=text.index('{',start)+1; depth=1
    while depth:
        if text[end]=='{':depth+=1
        if text[end]=='}':depth-=1
        end+=1
    return text[start:end]
helpers=[]
for name in ['PB32WideLess','PB32WideSubtract']:
    body=function(source,name);helpers.append(body);source=source.replace(body+'\n','',1)
needle='SIRIUS_PB32_INLINE uint PB32Divide(uint a, uint b) {'
assert source.count(needle)==1
source=source.replace(needle,'\n'.join(helpers)+'\n\n'+(W/'checked_estimate.h').read_text()+'\n'+needle)
needle='    // Normalization gives d <= remainder < 2*d.'
assert source.count(needle)==1
source=source.replace(needle,'''    double estimate = (double(remainder) * double(67108864)) / double(bb.significand);
    PB32Estimate checked = PB32CheckedDivide(remainder, bb.significand, estimate);
    if (checked.valid != 0u) return PB32Round(sign, exponent, checked.quotient | checked.sticky);
'''+needle)
needle='    uint remainder = 0u, root = 0u;'
assert source.count(needle)==1
source=source.replace(needle,'''    double scaled = double(aa.significand) * double(1u << shift);
#if defined(__cplusplus)
    double estimate = std::sqrt(scaled);
#else
    double estimate = sqrt(scaled);
#endif
    PB32Estimate checked = PB32CheckedSqrt(radicand, estimate);
    if (checked.valid != 0u) return PB32Round(0u, exponent, checked.quotient | checked.sticky);
'''+needle)
(W/'portable_binary32.h').write_text(source)
print(json.dumps({'private_header':str(W/'portable_binary32.h'),'integer_fallback_bodies_unchanged':True,'public_files_changed':False}))
