#pragma clang diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wmissing-prototypes"
#pragma clang diagnostic ignored "-Wmissing-braces"

#include <metal_stdlib>
#include <simd/simd.h>
#include <metal_atomic>

using namespace metal;

template<typename T, size_t Num>
struct spvUnsafeArray
{
    T elements[Num ? Num : 1];
    
    thread T& operator [] (size_t pos) thread
    {
        return elements[pos];
    }
    constexpr const thread T& operator [] (size_t pos) const thread
    {
        return elements[pos];
    }
    
    device T& operator [] (size_t pos) device
    {
        return elements[pos];
    }
    constexpr const device T& operator [] (size_t pos) const device
    {
        return elements[pos];
    }
    
    constexpr const constant T& operator [] (size_t pos) const constant
    {
        return elements[pos];
    }
    
    threadgroup T& operator [] (size_t pos) threadgroup
    {
        return elements[pos];
    }
    constexpr const threadgroup T& operator [] (size_t pos) const threadgroup
    {
        return elements[pos];
    }
};

// Implementation of the unsigned GLSL findMSB() function
template<typename T>
inline T spvFindUMSB(T x)
{
    return select(clz(T(0)) - (clz(x) + T(1)), T(-1), x == T(0));
}

template<typename T, typename U, typename V>
[[clang::optnone]] T spvMulExtended(V l, V r)
{
    return T{U(l * r), U(mulhi(l, r))};
}

struct RetainedTriple
{
    uint hi;
    uint lo;
    uint tail;
    uint error;
    uint valid;
};

struct StructuredBuffer
{
    uint _member0[1];
};

struct RWStructuredBuffer
{
    uint _member0[1];
};

struct PB32Finite
{
    uint significand;
    int exponent;
};

struct PB32Wide
{
    uint lo;
    uint hi;
};

struct PB32Product
{
    uint high;
    uint low;
    uint valid;
};

struct spvDescriptorSetBuffer0
{
    constant uint* spvBufferSizeConstants [[id(0)]];
    device StructuredBuffer* inputs [[id(1)]];
    device RWStructuredBuffer* outputs [[id(2)]];
};

static inline __attribute__((always_inline))
uint RPFromBits(uint bits)
{
    return bits;
}

static inline __attribute__((always_inline))
uint RPBits(uint value)
{
    return value;
}

static inline __attribute__((always_inline))
bool RPBounded(uint x, uint limit)
{
    return (RPBits(x) & 2147483647u) <= limit;
}

static inline __attribute__((always_inline))
uint RPZero()
{
    return RPFromBits(0u);
}

static inline __attribute__((always_inline))
RetainedTriple RTInvalid()
{
    RetainedTriple r;
    uint _385 = RPZero();
    r.hi = _385;
    r.lo = _385;
    r.tail = _385;
    r.error = _385;
    r.valid = 0u;
    return r;
}

static inline __attribute__((always_inline))
uint PB32Abs(uint a)
{
    return a & 2147483647u;
}

static inline __attribute__((always_inline))
bool PB32Nan(uint a)
{
    return PB32Abs(a) > 2139095040u;
}

static inline __attribute__((always_inline))
bool PB32Infinity(uint a)
{
    return PB32Abs(a) == 2139095040u;
}

static inline __attribute__((always_inline))
bool PB32Zero(uint a)
{
    return PB32Abs(a) == 0u;
}

static inline __attribute__((always_inline))
int PB32LeadingBit(uint a)
{
    return int(spvFindUMSB(a));
}

static inline __attribute__((always_inline))
PB32Finite PB32Unpack(uint a)
{
    uint e = (a >> uint(23)) & 255u;
    PB32Finite r;
    r.significand = a & 8388607u;
    r.exponent = int(e) - 127;
    if (e != 0u)
    {
        r.significand |= 8388608u;
    }
    else
    {
        r.exponent = -126;
        if (r.significand != 0u)
        {
            int _736 = PB32LeadingBit(r.significand);
            uint distance = uint(23 - _736);
            r.significand = r.significand << distance;
            r.exponent -= int(distance);
        }
    }
    return r;
}

static inline __attribute__((always_inline))
uint PB32ShiftJam(uint a, uint _distance)
{
    if (_distance == 0u)
    {
        return a;
    }
    if (_distance >= 32u)
    {
        uint _5119;
        if (a != 0u)
        {
            _5119 = 1u;
        }
        else
        {
            _5119 = 0u;
        }
        return _5119;
    }
    uint _5118;
    if ((a << (32u - _distance)) != 0u)
    {
        _5118 = 1u;
    }
    else
    {
        _5118 = 0u;
    }
    return (a >> _distance) | _5118;
}

static inline __attribute__((always_inline))
uint PB32Round(uint _sign, int exponent, uint extended)
{
    if (extended == 0u)
    {
        return _sign << uint(31);
    }
    uint _5120;
    int _5126;
    if (exponent < (-126))
    {
        _5126 = -126;
        _5120 = PB32ShiftJam(extended, uint((-126) - exponent));
    }
    else
    {
        _5126 = exponent;
        _5120 = extended;
    }
    uint residue = _5120 & 7u;
    uint sig = _5120 >> uint(3);
    bool _5121;
    if (residue > 4u)
    {
        _5121 = true;
    }
    else
    {
        bool _5122;
        if (residue == 4u)
        {
            _5122 = (sig & 1u) != 0u;
        }
        else
        {
            _5122 = false;
        }
        _5121 = _5122;
    }
    uint _5123;
    if (_5121)
    {
        _5123 = sig + 1u;
    }
    else
    {
        _5123 = sig;
    }
    int _5128;
    uint _5129;
    if (_5123 >= 16777216u)
    {
        _5129 = _5123 >> uint(1);
        _5128 = _5126 + 1;
    }
    else
    {
        _5129 = _5123;
        _5128 = _5126;
    }
    if (_5128 > 127)
    {
        return (_sign << uint(31)) | 2139095040u;
    }
    if (_5129 == 0u)
    {
        return _sign << uint(31);
    }
    bool _5130;
    if (_5128 == (-126))
    {
        _5130 = _5129 < 8388608u;
    }
    else
    {
        _5130 = false;
    }
    uint _5132;
    if (_5130)
    {
        _5132 = 0u;
    }
    else
    {
        _5132 = uint(_5128 + 127);
    }
    return ((_sign << uint(31)) | (_5132 << uint(23))) | (_5129 & 8388607u);
}

static inline __attribute__((always_inline))
uint PB32Add(uint a, uint b)
{
    bool _5106;
    if (PB32Nan(a))
    {
        _5106 = true;
    }
    else
    {
        _5106 = PB32Nan(b);
    }
    if (_5106)
    {
        return 2143289344u;
    }
    bool _621 = PB32Infinity(a);
    bool _5107;
    if (_621)
    {
        _5107 = true;
    }
    else
    {
        _5107 = PB32Infinity(b);
    }
    if (_5107)
    {
        bool _5115;
        if (_621)
        {
            _5115 = PB32Infinity(b);
        }
        else
        {
            _5115 = false;
        }
        bool _5116;
        if (_5115)
        {
            _5116 = ((a ^ b) >> uint(31)) != 0u;
        }
        else
        {
            _5116 = false;
        }
        if (_5116)
        {
            return 2143289344u;
        }
        uint _5117;
        if (_621)
        {
            _5117 = a;
        }
        else
        {
            _5117 = b;
        }
        return _5117;
    }
    bool _663 = PB32Zero(a);
    bool _5108;
    if (_663)
    {
        _5108 = PB32Zero(b);
    }
    else
    {
        _5108 = false;
    }
    if (_5108)
    {
        return (a & b) & 2147483648u;
    }
    if (_663)
    {
        return b;
    }
    if (PB32Zero(b))
    {
        return a;
    }
    uint _5109;
    uint _5110;
    if (PB32Abs(a) < PB32Abs(b))
    {
        _5110 = a;
        _5109 = b;
    }
    else
    {
        _5110 = b;
        _5109 = a;
    }
    PB32Finite aa = PB32Unpack(_5109);
    PB32Finite bb = PB32Unpack(_5110);
    uint ax = aa.significand << uint(3);
    int _764 = aa.exponent;
    uint bx = PB32ShiftJam(bb.significand << uint(3), uint(_764 - bb.exponent));
    int _5111;
    uint _5113;
    if (((_5109 ^ _5110) >> uint(31)) == 0u)
    {
        uint result = ax + bx;
        int _5112;
        uint _5114;
        if ((result & 134217728u) != 0u)
        {
            _5114 = PB32ShiftJam(result, 1u);
            _5112 = _764 + 1;
        }
        else
        {
            _5114 = result;
            _5112 = _764;
        }
        _5113 = _5114;
        _5111 = _5112;
    }
    else
    {
        uint result_1 = ax - bx;
        if (result_1 == 0u)
        {
            return 0u;
        }
        uint distance = uint(26 - PB32LeadingBit(result_1));
        _5113 = result_1 << distance;
        _5111 = _764 - int(distance);
    }
    return PB32Round(_5109 >> uint(31), _5111, _5113);
}

static inline __attribute__((always_inline))
uint OrderedAdd32(uint a, uint b)
{
    return PB32Add(a, b);
}

static inline __attribute__((always_inline))
uint PB32Negate(uint a)
{
    return a ^ 2147483648u;
}

static inline __attribute__((always_inline))
uint RPScalarNegate(uint value)
{
    return PB32Negate(value);
}

static inline __attribute__((always_inline))
void RetainedSum32(uint a, uint b, thread uint& hi, thread uint& lo)
{
    uint _5089;
    uint _5090;
    if ((RPBits(a) & 2147483647u) < (RPBits(b) & 2147483647u))
    {
        _5090 = a;
        _5089 = b;
    }
    else
    {
        _5090 = b;
        _5089 = a;
    }
    uint i;
    uint _5100;
    uint _5101;
    uint _5091 = 0u;
    uint _5092 = _5089;
    uint _5093 = _5090;
    for (;;)
    {
        bool _504_ladder_break = false;
        do
        {
            if (!(_5091 < 3u))
            {
                _504_ladder_break = true;
                break;
            }
            uint sum = OrderedAdd32(_5092, _5093);
            if (_5091 == 0u)
            {
                hi = sum;
                _5101 = RPScalarNegate(_5089);
                _5100 = sum;
            }
            else
            {
                uint _5095;
                uint _5096;
                if (_5091 == 1u)
                {
                    _5096 = RPScalarNegate(sum);
                    _5095 = _5090;
                }
                else
                {
                    lo = sum;
                    _5096 = _5093;
                    _5095 = _5092;
                }
                _5101 = _5096;
                _5100 = _5095;
            }
            break;
        } while(false);
        if (_504_ladder_break)
        {
            break;
        }
        i = _5091 + 1u;
        _5093 = _5101;
        _5092 = _5100;
        _5091 = i;
        continue;
    }
}

static inline __attribute__((always_inline))
RetainedTriple RTFinishAccumulation(uint hi, uint lo, uint tail, uint error)
{
    bool _5063;
    if (!RPBounded(error, 2071986176u))
    {
        _5063 = true;
    }
    else
    {
        _5063 = (RPBits(error) & 2147483648u) != 0u;
    }
    if (_5063)
    {
        return RTInvalid();
    }
    RetainedTriple _471 = RTInvalid();
    RetainedTriple r = _471;
    uint _473 = RPZero();
    uint i;
    uint upper;
    uint lower;
    uint _5077;
    uint _5081;
    uint _5085;
    uint _5088;
    uint _5064 = 0u;
    uint _5065 = lo;
    uint _5066 = tail;
    uint _5067 = _473;
    uint _5069 = _473;
    for (;;)
    {
        bool _434_ladder_break = false;
        do
        {
            if (!(_5064 < 5u))
            {
                _434_ladder_break = true;
                break;
            }
            RetainedSum32(_5065, _5066, upper, lower);
            if (_5064 == 0u)
            {
                _5088 = _5069;
                _5085 = lower;
                _5081 = upper;
                _5077 = hi;
            }
            else
            {
                uint _5072;
                uint _5078;
                uint _5082;
                uint _5086;
                if (_5064 == 1u)
                {
                    _5086 = _5067;
                    _5082 = _5067;
                    _5078 = lower;
                    _5072 = upper;
                }
                else
                {
                    uint _5079;
                    uint _5083;
                    uint _5087;
                    if (_5064 == 2u)
                    {
                        _5087 = lower;
                        _5083 = upper;
                        _5079 = _5069;
                    }
                    else
                    {
                        uint _5080;
                        uint _5084;
                        if (_5064 == 3u)
                        {
                            r.hi = upper;
                            _5084 = _5067;
                            _5080 = lower;
                        }
                        else
                        {
                            r.lo = upper;
                            r.tail = lower;
                            _5084 = _5066;
                            _5080 = _5065;
                        }
                        _5087 = _5067;
                        _5083 = _5084;
                        _5079 = _5080;
                    }
                    _5086 = _5087;
                    _5082 = _5083;
                    _5078 = _5079;
                    _5072 = _5069;
                }
                _5088 = _5072;
                _5085 = _5086;
                _5081 = _5082;
                _5077 = _5078;
            }
            break;
        } while(false);
        if (_434_ladder_break)
        {
            break;
        }
        i = _5064 + 1u;
        _5069 = _5088;
        _5067 = _5085;
        _5066 = _5081;
        _5065 = _5077;
        _5064 = i;
        continue;
    }
    r.error = error;
    r.valid = 1u;
    if (!RPBounded(r.hi, 2071986176u))
    {
        return _471;
    }
    return r;
}

static inline __attribute__((always_inline))
RetainedTriple RTPack(uint hi, uint lo, uint tail, uint error)
{
    bool _5061;
    if (!RPBounded(hi, 2071986176u))
    {
        _5061 = true;
    }
    else
    {
        _5061 = !RPBounded(lo, 2071986176u);
    }
    bool _5062;
    if (_5061)
    {
        _5062 = true;
    }
    else
    {
        _5062 = !RPBounded(tail, 2071986176u);
    }
    if (_5062)
    {
        return RTInvalid();
    }
    return RTFinishAccumulation(hi, lo, tail, error);
}

static inline __attribute__((always_inline))
RetainedTriple ReadOriginal(uint index, uint row, device StructuredBuffer& inputs, constant uint& inputsBufferSize)
{
    uint base = (1u + (225u * row)) + (5u * index);
    RetainedTriple value = RTPack(RPFromBits(inputs._member0[base]), RPFromBits(inputs._member0[base + 1u]), RPFromBits(inputs._member0[base + 2u]), RPFromBits(inputs._member0[base + 3u]));
    if (inputs._member0[base + 4u] != 1u)
    {
        value.valid = 0u;
    }
    return value;
}

static inline __attribute__((always_inline))
uint RPOne()
{
    return RPFromBits(1065353216u);
}

static inline __attribute__((always_inline))
bool PB32Equal(uint a, uint b)
{
    bool _5139;
    if (!PB32Nan(a))
    {
        _5139 = !PB32Nan(b);
    }
    else
    {
        _5139 = false;
    }
    bool _5140;
    if (_5139)
    {
        bool _5141;
        if (a == b)
        {
            _5141 = true;
        }
        else
        {
            bool _5142;
            if (PB32Zero(a))
            {
                _5142 = PB32Zero(b);
            }
            else
            {
                _5142 = false;
            }
            _5141 = _5142;
        }
        _5140 = _5141;
    }
    else
    {
        _5140 = false;
    }
    return _5140;
}

static inline __attribute__((always_inline))
bool RPScalarEqual(uint a, uint b)
{
    return PB32Equal(a, b);
}

static inline __attribute__((always_inline))
bool RTExact(thread const RetainedTriple& a, uint x)
{
    bool _5135;
    if (a.valid != 0u)
    {
        _5135 = RPScalarEqual(a.hi, x);
    }
    else
    {
        _5135 = false;
    }
    bool _5136;
    if (_5135)
    {
        _5136 = RPScalarEqual(a.lo, RPZero());
    }
    else
    {
        _5136 = false;
    }
    bool _5137;
    if (_5136)
    {
        _5137 = RPScalarEqual(a.tail, RPZero());
    }
    else
    {
        _5137 = false;
    }
    bool _5138;
    if (_5137)
    {
        _5138 = RPScalarEqual(a.error, RPZero());
    }
    else
    {
        _5138 = false;
    }
    return _5138;
}

static inline __attribute__((always_inline))
uint RPAbs(uint x)
{
    return RPFromBits(RPBits(x) & 2147483647u);
}

static inline __attribute__((always_inline))
uint RPFailureBound()
{
    return RPFromBits(2130706432u);
}

static inline __attribute__((always_inline))
uint RPUpAdd(uint a, uint b)
{
    bool _5143;
    if (!RPBounded(a, 2105540608u))
    {
        _5143 = true;
    }
    else
    {
        _5143 = !RPBounded(b, 2105540608u);
    }
    if (_5143)
    {
        return RPFailureBound();
    }
    uint _1324 = RPZero();
    if (RPScalarEqual(a, _1324))
    {
        return b;
    }
    if (RPScalarEqual(b, _1324))
    {
        return a;
    }
    return RPFromBits(RPBits(OrderedAdd32(a, b)) + 1u);
}

static inline __attribute__((always_inline))
uint RTCenterUpper(thread const RetainedTriple& a)
{
    return RPUpAdd(RPUpAdd(RPAbs(a.hi), RPAbs(a.lo)), RPAbs(a.tail));
}

static inline __attribute__((always_inline))
uint RTAbsoluteUpper(thread const RetainedTriple& a)
{
    return RPUpAdd(RTCenterUpper(a), a.error);
}

static inline __attribute__((always_inline))
RetainedTriple RTNegate(thread const RetainedTriple& a)
{
    RetainedTriple _1412 = a;
    _1412.hi = RPFromBits(RPBits(_1412.hi) ^ 2147483648u);
    _1412.lo = RPFromBits(RPBits(_1412.lo) ^ 2147483648u);
    _1412.tail = RPFromBits(RPBits(_1412.tail) ^ 2147483648u);
    return _1412;
}

static inline __attribute__((always_inline))
void WriteStored(uint base, thread const RetainedTriple& r, device RWStructuredBuffer& outputs, constant uint& outputsBufferSize)
{
    outputs._member0[base] = RPBits(r.hi);
    outputs._member0[base + 1u] = RPBits(r.lo);
    outputs._member0[base + 2u] = RPBits(r.tail);
    outputs._member0[base + 3u] = RPBits(r.error);
    outputs._member0[base + 4u] = r.valid;
}

static inline __attribute__((always_inline))
void RejectRetainedProgram(threadgroup uint& retainedStatus)
{
    uint _1665 = atomic_fetch_or_explicit((threadgroup atomic_uint*)&retainedStatus, 1u, memory_order_relaxed);
}

static inline __attribute__((always_inline))
RetainedTriple RTConstant(uint x)
{
    if (!RPBounded(x, 2071986176u))
    {
        return RTInvalid();
    }
    RetainedTriple r;
    uint _5184;
    if ((RPBits(x) & 2147483647u) == 0u)
    {
        _5184 = RPZero();
    }
    else
    {
        _5184 = x;
    }
    r.hi = _5184;
    uint _1828 = RPZero();
    r.lo = _1828;
    r.tail = _1828;
    r.error = _1828;
    r.valid = 1u;
    return r;
}

static inline __attribute__((always_inline))
uint ScratchAddress(uint index, uint component, uint rowBase, uint prefix)
{
    return (index * 5u) + component;
}

static inline __attribute__((always_inline))
uint ScratchRead(uint index, uint component, uint rowBase, uint prefix, threadgroup spvUnsafeArray<uint, 1810>& retainedScratch)
{
    return retainedScratch[ScratchAddress(index, component, rowBase, prefix)];
}

static inline __attribute__((always_inline))
RetainedTriple ReadValue(uint index, uint rowBase, threadgroup spvUnsafeArray<uint, 1810>& retainedScratch)
{
    RetainedTriple r;
    r.hi = RPFromBits(ScratchRead(index, 0u, rowBase, 204u, retainedScratch));
    r.lo = RPFromBits(ScratchRead(index, 1u, rowBase, 204u, retainedScratch));
    r.tail = RPFromBits(ScratchRead(index, 2u, rowBase, 204u, retainedScratch));
    r.error = RPFromBits(ScratchRead(index, 3u, rowBase, 204u, retainedScratch));
    r.valid = ScratchRead(index, 4u, rowBase, 204u, retainedScratch);
    return r;
}

static inline __attribute__((always_inline))
bool PB32Less(uint a, uint b)
{
    bool _5316;
    if (PB32Nan(a))
    {
        _5316 = true;
    }
    else
    {
        _5316 = PB32Nan(b);
    }
    bool _5317;
    if (_5316)
    {
        _5317 = true;
    }
    else
    {
        _5317 = PB32Equal(a, b);
    }
    if (_5317)
    {
        return false;
    }
    uint sa = a >> uint(31);
    bool _5318;
    if (sa != (b >> uint(31)))
    {
        _5318 = sa != 0u;
    }
    else
    {
        bool _5319;
        if (sa != 0u)
        {
            _5319 = PB32Abs(a) > PB32Abs(b);
        }
        else
        {
            _5319 = PB32Abs(a) < PB32Abs(b);
        }
        _5318 = _5319;
    }
    return _5318;
}

static inline __attribute__((always_inline))
bool PB32LessEqual(uint a, uint b)
{
    bool _5315;
    if (PB32Equal(a, b))
    {
        _5315 = true;
    }
    else
    {
        _5315 = PB32Less(a, b);
    }
    return _5315;
}

static inline __attribute__((always_inline))
bool RPScalarLessEqual(uint a, uint b)
{
    return PB32LessEqual(a, b);
}

static inline __attribute__((always_inline))
void RetainedSum32_1(uint a, uint b, thread uint& hi, thread uint& lo)
{
    uint _5323;
    uint _5324;
    if ((RPBits(a) & 2147483647u) < (RPBits(b) & 2147483647u))
    {
        _5324 = a;
        _5323 = b;
    }
    else
    {
        _5324 = b;
        _5323 = a;
    }
    uint i;
    uint _5334;
    uint _5335;
    uint _5325 = 0u;
    uint _5326 = _5323;
    uint _5327 = _5324;
    for (;;)
    {
        bool _2292_ladder_break = false;
        do
        {
            if (!(_5325 < 3u))
            {
                _2292_ladder_break = true;
                break;
            }
            uint sum = OrderedAdd32(_5326, _5327);
            if (_5325 == 0u)
            {
                hi = sum;
                _5335 = RPScalarNegate(_5323);
                _5334 = sum;
            }
            else
            {
                uint _5329;
                uint _5330;
                if (_5325 == 1u)
                {
                    _5330 = RPScalarNegate(sum);
                    _5329 = _5324;
                }
                else
                {
                    lo = sum;
                    _5330 = _5327;
                    _5329 = _5326;
                }
                _5335 = _5330;
                _5334 = _5329;
            }
            break;
        } while(false);
        if (_2292_ladder_break)
        {
            break;
        }
        i = _5325 + 1u;
        _5327 = _5335;
        _5326 = _5334;
        _5325 = i;
        continue;
    }
}

static inline __attribute__((always_inline))
bool RPScalarGreaterEqual(uint a, uint b)
{
    return RPScalarLessEqual(b, a);
}

static inline __attribute__((always_inline))
uint RTAbsoluteLower(thread const RetainedTriple& a)
{
    uint tail = RPUpAdd(RPUpAdd(RPAbs(a.lo), RPAbs(a.tail)), a.error);
    bool _5320;
    if (a.valid == 0u)
    {
        _5320 = true;
    }
    else
    {
        _5320 = RPScalarLessEqual(RPAbs(a.hi), tail);
    }
    if (_5320)
    {
        return RPZero();
    }
    uint difference;
    uint correction;
    RetainedSum32_1(RPAbs(a.hi), RPScalarNegate(tail), difference, correction);
    uint bits = RPBits(difference);
    uint _2368 = RPZero();
    uint _5321;
    if (RPScalarGreaterEqual(correction, _2368))
    {
        _5321 = difference;
    }
    else
    {
        uint _5322;
        if (bits == 0u)
        {
            _5322 = _2368;
        }
        else
        {
            _5322 = RPFromBits(bits - 1u);
        }
        _5321 = _5322;
    }
    return _5321;
}

static inline __attribute__((always_inline))
int RPFloorExponent(uint positive)
{
    uint bits = RPBits(positive);
    uint exponent = bits >> uint(23);
    int _5341;
    if (exponent != 0u)
    {
        _5341 = int(exponent) - 127;
    }
    else
    {
        _5341 = int(spvFindUMSB(bits)) - 149;
    }
    return _5341;
}

static inline __attribute__((always_inline))
uint RPMantissa(uint positive)
{
    uint bits = RPBits(positive);
    uint _5342;
    if ((bits >> uint(23)) != 0u)
    {
        _5342 = bits & 8388607u;
    }
    else
    {
        _5342 = (bits << uint(23 - int(spvFindUMSB(bits)))) & 8388607u;
    }
    return RPFromBits(1065353216u | _5342);
}

static inline __attribute__((always_inline))
uint PB32WideShiftJam(thread const PB32Wide& a, uint _distance)
{
    if (_distance == 0u)
    {
        return a.lo;
    }
    if (_distance < 32u)
    {
        uint _2612 = 32u - _distance;
        uint _5353;
        if ((a.lo << _2612) != 0u)
        {
            _5353 = 1u;
        }
        else
        {
            _5353 = 0u;
        }
        return ((a.lo >> _distance) | (a.hi << _2612)) | _5353;
    }
    if (_distance == 32u)
    {
        uint _5352;
        if (a.lo != 0u)
        {
            _5352 = 1u;
        }
        else
        {
            _5352 = 0u;
        }
        return a.hi | _5352;
    }
    if (_distance < 64u)
    {
        uint _5351;
        if (((a.hi << (64u - _distance)) | a.lo) != 0u)
        {
            _5351 = 1u;
        }
        else
        {
            _5351 = 0u;
        }
        return (a.hi >> (_distance - 32u)) | _5351;
    }
    uint _5350;
    if ((a.lo | a.hi) != 0u)
    {
        _5350 = 1u;
    }
    else
    {
        _5350 = 0u;
    }
    return _5350;
}

static inline __attribute__((always_inline))
uint PB32Multiply(uint a, uint b)
{
    uint sign = (a ^ b) >> uint(31);
    bool _5343;
    if (PB32Nan(a))
    {
        _5343 = true;
    }
    else
    {
        _5343 = PB32Nan(b);
    }
    if (_5343)
    {
        return 2143289344u;
    }
    bool _5344;
    if (PB32Infinity(a))
    {
        _5344 = true;
    }
    else
    {
        _5344 = PB32Infinity(b);
    }
    if (_5344)
    {
        bool _5348;
        if (PB32Zero(a))
        {
            _5348 = true;
        }
        else
        {
            _5348 = PB32Zero(b);
        }
        uint _5349;
        if (_5348)
        {
            _5349 = 2143289344u;
        }
        else
        {
            _5349 = (sign << uint(31)) | 2139095040u;
        }
        return _5349;
    }
    bool _5345;
    if (PB32Zero(a))
    {
        _5345 = true;
    }
    else
    {
        _5345 = PB32Zero(b);
    }
    if (_5345)
    {
        return sign << uint(31);
    }
    PB32Finite aa = PB32Unpack(a);
    PB32Finite bb = PB32Unpack(b);
    PB32Wide _2559 = spvMulExtended<PB32Wide, uint>(aa.significand, bb.significand);
    bool upper = (_2559.hi & 32768u) != 0u;
    uint _5346;
    if (upper)
    {
        _5346 = 21u;
    }
    else
    {
        _5346 = 20u;
    }
    PB32Wide _2477 = _2559;
    int _5347;
    if (upper)
    {
        _5347 = 1;
    }
    else
    {
        _5347 = 0;
    }
    return PB32Round(sign, (aa.exponent + bb.exponent) + _5347, PB32WideShiftJam(_2477, _5346));
}

static inline __attribute__((always_inline))
uint OrderedMultiply32(uint a, uint b)
{
    return PB32Multiply(a, b);
}

static inline __attribute__((always_inline))
uint PB32Sqrt(uint a)
{
    if (PB32Nan(a))
    {
        return 2143289344u;
    }
    if (PB32Zero(a))
    {
        return a;
    }
    if ((a >> uint(31)) != 0u)
    {
        return 2143289344u;
    }
    if (PB32Infinity(a))
    {
        return a;
    }
    PB32Finite aa = PB32Unpack(a);
    int _2744 = aa.exponent;
    uint odd = uint(_2744) & 1u;
    uint shift = 29u + odd;
    PB32Wide radicand;
    uint _2753 = aa.significand;
    radicand.lo = _2753 << shift;
    radicand.hi = _2753 >> (32u - shift);
    uint _5355;
    uint _5356;
    _5356 = 0u;
    _5355 = 0u;
    int i;
    uint _5363;
    uint _5364;
    int _5354 = 26;
    for (;;)
    {
        bool _2714_ladder_break = false;
        do
        {
            if (!(_5354 >= 0))
            {
                _2714_ladder_break = true;
                break;
            }
            uint bit = uint(_5354) * 2u;
            uint _5359;
            if (bit >= 32u)
            {
                _5359 = (radicand.hi >> (bit - 32u)) & 3u;
            }
            else
            {
                _5359 = (radicand.lo >> bit) & 3u;
            }
            uint remainder = (_5355 << uint(2)) | _5359;
            uint trial = (_5356 << uint(2)) | 1u;
            uint root = _5356 << uint(1);
            if (remainder >= trial)
            {
                _5364 = root | 1u;
                _5363 = remainder - trial;
            }
            else
            {
                _5364 = root;
                _5363 = remainder;
            }
            break;
        } while(false);
        if (_2714_ladder_break)
        {
            break;
        }
        i = _5354 - 1;
        _5356 = _5364;
        _5355 = _5363;
        _5354 = i;
        continue;
    }
    uint _5357;
    if (_5355 != 0u)
    {
        _5357 = _5356 | 1u;
    }
    else
    {
        _5357 = _5356;
    }
    return PB32Round(0u, (_2744 - int(odd)) / 2, _5357);
}

static inline __attribute__((always_inline))
uint OrderedSqrt32(uint value)
{
    return PB32Sqrt(value);
}

static inline __attribute__((always_inline))
uint RPPowerOfTwo(int e)
{
    if (e > 120)
    {
        return RPFailureBound();
    }
    if (e < (-149))
    {
        return RPFromBits(1u);
    }
    uint _5365;
    if (e >= (-126))
    {
        _5365 = RPFromBits(uint(e + 127) << uint(23));
    }
    else
    {
        _5365 = RPFromBits(1u << uint(e + 149));
    }
    return _5365;
}

static inline __attribute__((always_inline))
uint RPSqrtSeed(uint positive)
{
    int exponent = RPFloorExponent(positive);
    int odd = exponent & 1;
    uint mantissa = RPMantissa(positive);
    uint _5340;
    if (odd != 0)
    {
        _5340 = OrderedMultiply32(mantissa, RPFromBits(1073741824u));
    }
    else
    {
        _5340 = mantissa;
    }
    uint seed = OrderedSqrt32(_5340);
    if (!RPBounded(seed, 1082130432u))
    {
        return RPZero();
    }
    return OrderedMultiply32(seed, RPPowerOfTwo((exponent - odd) / 2));
}

static inline __attribute__((always_inline))
uint RPRatioUpper(uint numerator, uint denominator)
{
    bool _5366;
    if (!RPBounded(numerator, 2071986176u))
    {
        _5366 = true;
    }
    else
    {
        _5366 = !RPBounded(denominator, 2071986176u);
    }
    bool _5367;
    if (_5366)
    {
        _5367 = true;
    }
    else
    {
        _5367 = RPScalarLessEqual(denominator, RPZero());
    }
    if (_5367)
    {
        return RPFailureBound();
    }
    uint _2929 = RPZero();
    if (RPScalarEqual(numerator, _2929))
    {
        return _2929;
    }
    return RPPowerOfTwo((RPFloorExponent(numerator) - RPFloorExponent(denominator)) + 1);
}

static inline __attribute__((always_inline))
bool RPScalarLess(uint a, uint b)
{
    return PB32Less(a, b);
}

static inline __attribute__((always_inline))
bool RPScalarGreater(uint a, uint b)
{
    return RPScalarLess(b, a);
}

static inline __attribute__((always_inline))
uint PB32Divide(uint a, uint b)
{
    uint sign = (a ^ b) >> uint(31);
    bool _5371;
    if (PB32Nan(a))
    {
        _5371 = true;
    }
    else
    {
        _5371 = PB32Nan(b);
    }
    if (_5371)
    {
        return 2143289344u;
    }
    bool _3082 = PB32Infinity(a);
    bool _5372;
    if (_3082)
    {
        _5372 = PB32Infinity(b);
    }
    else
    {
        _5372 = false;
    }
    bool _5373;
    if (_5372)
    {
        _5373 = true;
    }
    else
    {
        bool _5374;
        if (PB32Zero(a))
        {
            _5374 = PB32Zero(b);
        }
        else
        {
            _5374 = false;
        }
        _5373 = _5374;
    }
    if (_5373)
    {
        return 2143289344u;
    }
    bool _5375;
    if (_3082)
    {
        _5375 = true;
    }
    else
    {
        _5375 = PB32Zero(b);
    }
    if (_5375)
    {
        return (sign << uint(31)) | 2139095040u;
    }
    bool _5376;
    if (PB32Zero(a))
    {
        _5376 = true;
    }
    else
    {
        _5376 = PB32Infinity(b);
    }
    if (_5376)
    {
        return sign << uint(31);
    }
    PB32Finite aa = PB32Unpack(a);
    PB32Finite bb = PB32Unpack(b);
    int exponent = aa.exponent - bb.exponent;
    uint remainder = aa.significand;
    uint _3132 = bb.significand;
    uint _5377;
    int _5383;
    if (remainder < _3132)
    {
        _5383 = exponent - 1;
        _5377 = remainder << uint(1);
    }
    else
    {
        _5383 = exponent;
        _5377 = remainder;
    }
    uint _5379;
    uint _5380;
    _5380 = 1u;
    _5379 = _5377 - _3132;
    uint _3163;
    uint _3167;
    uint chunk;
    uint _5378 = 0u;
    for (;;)
    {
        bool _3057_ladder_break = false;
        do
        {
            if (!(_5378 < 4u))
            {
                _3057_ladder_break = true;
                break;
            }
            uint _5386;
            if (_5378 == 3u)
            {
                _5386 = 2u;
            }
            else
            {
                _5386 = 8u;
            }
            uint numerator = _5379 << _5386;
            uint digit = numerator / _3132;
            _3163 = numerator - (digit * _3132);
            _3167 = (_5380 << _5386) | digit;
            break;
        } while(false);
        if (_3057_ladder_break)
        {
            break;
        }
        chunk = _5378 + 1u;
        _5380 = _3167;
        _5379 = _3163;
        _5378 = chunk;
        continue;
    }
    uint _5384;
    if (_5379 != 0u)
    {
        _5384 = _5380 | 1u;
    }
    else
    {
        _5384 = _5380;
    }
    return PB32Round(sign, _5383, _5384);
}

static inline __attribute__((always_inline))
uint OrderedDivide32(uint a, uint b)
{
    return PB32Divide(a, b);
}

static inline __attribute__((always_inline))
uint RPDivideSeed(uint a, uint b)
{
    uint _2988 = RPZero();
    if (RPScalarEqual(a, _2988))
    {
        return _2988;
    }
    uint aa = RPAbs(a);
    uint bb = RPAbs(b);
    int exponent = RPFloorExponent(aa) - RPFloorExponent(bb);
    bool _5368;
    if (exponent < (-151))
    {
        _5368 = true;
    }
    else
    {
        _5368 = exponent > 120;
    }
    if (_5368)
    {
        return _2988;
    }
    uint seed = OrderedDivide32(RPMantissa(aa), RPMantissa(bb));
    if (!RPBounded(seed, 1082130432u))
    {
        return _2988;
    }
    int _5369;
    uint _5370;
    if (exponent < (-126))
    {
        _5370 = OrderedMultiply32(seed, RPFromBits(8388608u));
        _5369 = exponent + 126;
    }
    else
    {
        _5370 = seed;
        _5369 = exponent;
    }
    return RPFromBits(RPBits(OrderedMultiply32(_5370, RPPowerOfTwo(_5369))) ^ ((RPBits(a) ^ RPBits(b)) & 2147483648u));
}

static inline __attribute__((always_inline))
uint PB32FromUnsigned(uint a)
{
    if (a == 0u)
    {
        return 0u;
    }
    int exponent = PB32LeadingBit(a);
    uint _5418;
    if (exponent <= 26)
    {
        _5418 = a << uint(26 - exponent);
    }
    else
    {
        _5418 = PB32ShiftJam(a, uint(exponent - 26));
    }
    return PB32Round(0u, exponent, _5418);
}

static inline __attribute__((always_inline))
uint PB32FromSignedBits(uint a)
{
    uint sign = a >> uint(31);
    uint _5417;
    if (sign != 0u)
    {
        _5417 = 0u - a;
    }
    else
    {
        _5417 = a;
    }
    return PB32FromUnsigned(_5417) | (sign << uint(31));
}

static inline __attribute__((always_inline))
uint RPFromSigned(int value)
{
    return PB32FromSignedBits(uint(value));
}

static inline __attribute__((always_inline))
bool PB32WideLess(thread const PB32Wide& a, thread const PB32Wide& b)
{
    bool _5433;
    if (a.hi != b.hi)
    {
        _5433 = a.hi < b.hi;
    }
    else
    {
        _5433 = a.lo < b.lo;
    }
    return _5433;
}

static inline __attribute__((always_inline))
PB32Wide PB32WideSubtract(thread const PB32Wide& a, thread const PB32Wide& b)
{
    PB32Wide r;
    r.lo = a.lo - b.lo;
    uint _5434;
    if (a.lo < b.lo)
    {
        _5434 = 1u;
    }
    else
    {
        _5434 = 0u;
    }
    r.hi = (a.hi - b.hi) - _5434;
    return r;
}

static inline __attribute__((always_inline))
uint PB32WideShiftJam_1(thread const PB32Wide& a, uint _distance)
{
    if (_distance == 0u)
    {
        return a.lo;
    }
    if (_distance < 32u)
    {
        uint _3919 = 32u - _distance;
        uint _5441;
        if ((a.lo << _3919) != 0u)
        {
            _5441 = 1u;
        }
        else
        {
            _5441 = 0u;
        }
        return ((a.lo >> _distance) | (a.hi << _3919)) | _5441;
    }
    if (_distance == 32u)
    {
        uint _5440;
        if (a.lo != 0u)
        {
            _5440 = 1u;
        }
        else
        {
            _5440 = 0u;
        }
        return a.hi | _5440;
    }
    if (_distance < 64u)
    {
        uint _5439;
        if (((a.hi << (64u - _distance)) | a.lo) != 0u)
        {
            _5439 = 1u;
        }
        else
        {
            _5439 = 0u;
        }
        return (a.hi >> (_distance - 32u)) | _5439;
    }
    uint _5438;
    if ((a.lo | a.hi) != 0u)
    {
        _5438 = 1u;
    }
    else
    {
        _5438 = 0u;
    }
    return _5438;
}

static inline __attribute__((always_inline))
uint PB32RoundWide(uint _sign, int leastBitExponent, thread const PB32Wide& magnitude)
{
    if ((magnitude.lo | magnitude.hi) == 0u)
    {
        return _sign << uint(31);
    }
    int _5435;
    if (magnitude.hi != 0u)
    {
        _5435 = 32 + PB32LeadingBit(magnitude.hi);
    }
    else
    {
        _5435 = PB32LeadingBit(magnitude.lo);
    }
    uint _5437;
    if (_5435 > 26)
    {
        _5437 = PB32WideShiftJam_1(magnitude, uint(_5435 - 26));
    }
    else
    {
        _5437 = magnitude.lo << uint(26 - _5435);
    }
    return PB32Round(_sign, leastBitExponent + _5435, _5437);
}

static inline __attribute__((always_inline))
PB32Product PB32ProductResidual(uint a, uint b)
{
    PB32Product r;
    r.high = 0u;
    r.low = 0u;
    r.valid = 0u;
    bool _5427;
    if (PB32Abs(a) > 1568669696u)
    {
        _5427 = true;
    }
    else
    {
        _5427 = PB32Abs(b) > 1568669696u;
    }
    if (_5427)
    {
        return r;
    }
    r.high = PB32Multiply(a, b);
    bool _5428;
    if (PB32Zero(a))
    {
        _5428 = true;
    }
    else
    {
        _5428 = PB32Zero(b);
    }
    if (_5428)
    {
        r.valid = 1u;
        return r;
    }
    PB32Finite aa = PB32Unpack(a);
    PB32Finite bb = PB32Unpack(b);
    int sumExponent = aa.exponent + bb.exponent;
    PB32Wide _3705 = spvMulExtended<PB32Wide, uint>(aa.significand, bb.significand);
    PB32Wide rounded;
    rounded.lo = 0u;
    rounded.hi = 0u;
    if (!PB32Zero(r.high))
    {
        PB32Finite high = PB32Unpack(r.high);
        int shift = (high.exponent - sumExponent) + 23;
        bool _5429;
        if (shift < 23)
        {
            _5429 = true;
        }
        else
        {
            _5429 = shift > 25;
        }
        if (_5429)
        {
            r.high = 0u;
            return r;
        }
        uint _3732 = high.significand;
        rounded.lo = _3732 << uint(shift);
        rounded.hi = _3732 >> uint(32 - shift);
    }
    PB32Wide _3631 = _3705;
    PB32Wide _3632 = rounded;
    bool _3745 = PB32WideLess(_3631, _3632);
    PB32Wide _5430;
    if (_3745)
    {
        PB32Wide _3634 = rounded;
        PB32Wide _3635 = _3705;
        _5430 = PB32WideSubtract(_3634, _3635);
    }
    else
    {
        PB32Wide _3636 = _3705;
        PB32Wide _3637 = rounded;
        _5430 = PB32WideSubtract(_3636, _3637);
    }
    if ((_5430.lo | _5430.hi) != 0u)
    {
        uint _5431;
        if (_3745)
        {
            _5431 = 1u;
        }
        else
        {
            _5431 = 0u;
        }
        PB32Wide _3639 = _5430;
        r.low = PB32RoundWide(((a ^ b) >> uint(31)) ^ _5431, sumExponent - 46, _3639);
    }
    r.valid = 1u;
    return r;
}

static inline __attribute__((always_inline))
void RetainedProduct32(uint a, uint b, thread uint& hi, thread uint& lo)
{
    PB32Product product = PB32ProductResidual(RPBits(a), RPBits(b));
    bool _3999 = product.valid != 0u;
    uint _5425;
    if (_3999)
    {
        _5425 = product.high;
    }
    else
    {
        _5425 = 2139095040u;
    }
    hi = RPFromBits(_5425);
    uint _5426;
    if (_3999)
    {
        _5426 = product.low;
    }
    else
    {
        _5426 = 2139095040u;
    }
    lo = RPFromBits(_5426);
}

static inline __attribute__((always_inline))
RetainedTriple RTAccumulateTerm(thread const RetainedTriple& a, uint x)
{
    bool _5442;
    if (a.valid == 0u)
    {
        _5442 = true;
    }
    else
    {
        _5442 = !RPBounded(x, 2071986176u);
    }
    if (_5442)
    {
        return RTInvalid();
    }
    RetainedTriple _4102 = RTInvalid();
    RetainedTriple r = _4102;
    uint _5444;
    _5444 = x;
    uint _4148;
    uint i;
    uint upper;
    uint lower;
    uint _5443 = 0u;
    for (;;)
    {
        bool _4057_ladder_break = false;
        do
        {
            if (!(_5443 < 3u))
            {
                _4057_ladder_break = true;
                break;
            }
            bool _4110 = _5443 == 0u;
            uint _5450;
            if (_4110)
            {
                _5450 = a.hi;
            }
            else
            {
                uint _5451;
                if (_5443 == 1u)
                {
                    _5451 = a.lo;
                }
                else
                {
                    _5451 = a.tail;
                }
                _5450 = _5451;
            }
            RetainedSum32_1(_5450, _5444, upper, lower);
            if (_4110)
            {
                r.hi = upper;
            }
            else
            {
                if (_5443 == 1u)
                {
                    r.lo = upper;
                }
                else
                {
                    r.tail = upper;
                }
            }
            _4148 = lower;
            break;
        } while(false);
        if (_4057_ladder_break)
        {
            break;
        }
        i = _5443 + 1u;
        _5444 = _4148;
        _5443 = i;
        continue;
    }
    r.error = RPUpAdd(a.error, RPAbs(_5444));
    bool _5445;
    if (!RPBounded(r.hi, 2105540608u))
    {
        _5445 = true;
    }
    else
    {
        _5445 = !RPBounded(r.lo, 2105540608u);
    }
    bool _5446;
    if (_5445)
    {
        _5446 = true;
    }
    else
    {
        _5446 = !RPBounded(r.tail, 2105540608u);
    }
    bool _5447;
    if (_5446)
    {
        _5447 = true;
    }
    else
    {
        _5447 = !RPBounded(r.error, 2105540608u);
    }
    if (_5447)
    {
        return _4102;
    }
    r.valid = 1u;
    return r;
}

static inline __attribute__((always_inline))
RetainedTriple RTAccumulateProduct(thread const RetainedTriple& sum, uint a, uint b)
{
    bool _5419;
    if (sum.valid == 0u)
    {
        _5419 = true;
    }
    else
    {
        _5419 = !RPBounded(a, 1568669696u);
    }
    bool _5420;
    if (_5419)
    {
        _5420 = true;
    }
    else
    {
        _5420 = !RPBounded(b, 1568669696u);
    }
    if (_5420)
    {
        return RTInvalid();
    }
    uint _3592 = RPZero();
    bool _5421;
    if (RPScalarEqual(a, _3592))
    {
        _5421 = true;
    }
    else
    {
        _5421 = RPScalarEqual(b, _3592);
    }
    if (_5421)
    {
        return sum;
    }
    uint hi;
    uint lo;
    RetainedProduct32(a, b, hi, lo);
    RetainedTriple r = sum;
    uint i;
    uint _5422 = 0u;
    for (;;)
    {
        bool _3556_ladder_break = false;
        do
        {
            if (!(_5422 < 2u))
            {
                _3556_ladder_break = true;
                break;
            }
            uint _5423;
            if (_5422 == 0u)
            {
                _5423 = hi;
            }
            else
            {
                _5423 = lo;
            }
            RetainedTriple _3540 = r;
            r = RTAccumulateTerm(_3540, _5423);
            break;
        } while(false);
        if (_3556_ladder_break)
        {
            break;
        }
        i = _5422 + 1u;
        _5422 = i;
        continue;
    }
    if (r.valid != 0u)
    {
        r.error = RPUpAdd(r.error, RPFromBits(5u));
    }
    return r;
}

static inline __attribute__((always_inline))
uint RPUpMultiply(uint a, uint b)
{
    uint _4352 = RPZero();
    bool _5459;
    if (RPScalarEqual(a, _4352))
    {
        _5459 = true;
    }
    else
    {
        _5459 = RPScalarEqual(b, _4352);
    }
    if (_5459)
    {
        return _4352;
    }
    bool _5460;
    if (!RPBounded(a, 1585446912u))
    {
        _5460 = true;
    }
    else
    {
        _5460 = !RPBounded(b, 1585446912u);
    }
    if (_5460)
    {
        return RPFailureBound();
    }
    return RPFromBits(RPBits(OrderedMultiply32(a, b)) + 1u);
}

static inline __attribute__((always_inline))
RetainedTriple RTAccumulateRoundedProduct(thread const RetainedTriple& sum, uint a, uint b)
{
    bool _5456;
    if (sum.valid == 0u)
    {
        _5456 = true;
    }
    else
    {
        _5456 = !RPBounded(a, 1568669696u);
    }
    bool _5457;
    if (_5456)
    {
        _5457 = true;
    }
    else
    {
        _5457 = !RPBounded(b, 1568669696u);
    }
    if (_5457)
    {
        return RTInvalid();
    }
    uint _4320 = RPZero();
    bool _5458;
    if (RPScalarEqual(a, _4320))
    {
        _5458 = true;
    }
    else
    {
        _5458 = RPScalarEqual(b, _4320);
    }
    if (_5458)
    {
        return sum;
    }
    uint product = OrderedMultiply32(a, b);
    RetainedTriple r = RTAccumulateTerm(sum, product);
    if (r.valid != 0u)
    {
        r.error = RPUpAdd(r.error, RPUpAdd(RPUpMultiply(RPAbs(product), RPFromBits(872415232u)), RPFromBits(5u)));
    }
    return r;
}

static inline __attribute__((always_inline))
RetainedTriple RTMultiply(thread const RetainedTriple& a, thread const RetainedTriple& b)
{
    bool _5390;
    if (a.valid == 0u)
    {
        _5390 = true;
    }
    else
    {
        _5390 = b.valid == 0u;
    }
    if (_5390)
    {
        return RTInvalid();
    }
    uint _3410 = RPZero();
    bool _5391;
    if (RTExact(a, _3410))
    {
        _5391 = true;
    }
    else
    {
        _5391 = RTExact(b, _3410);
    }
    if (_5391)
    {
        return RTConstant(_3410);
    }
    uint _3422 = RPOne();
    if (RTExact(a, _3422))
    {
        return b;
    }
    if (RTExact(b, _3422))
    {
        return a;
    }
    uint _3431 = RPFromSigned(-1);
    if (RTExact(a, _3431))
    {
        return RTNegate(b);
    }
    if (RTExact(b, _3431))
    {
        return RTNegate(a);
    }
    RetainedTriple _5395;
    _5395 = RTConstant(_3410);
    RetainedTriple _3527;
    uint i;
    uint _5392 = 0u;
    for (;;)
    {
        bool _3357_ladder_break = false;
        do
        {
            if (!(_5392 < 3u))
            {
                _3357_ladder_break = true;
                break;
            }
            uint _5414;
            if (_5392 == 2u)
            {
                _5414 = a.lo;
            }
            else
            {
                _5414 = a.hi;
            }
            uint _5415;
            if (_5392 == 1u)
            {
                _5415 = b.lo;
            }
            else
            {
                _5415 = b.hi;
            }
            RetainedTriple _3335 = _5395;
            _3527 = RTAccumulateProduct(_3335, _5414, _5415);
            break;
        } while(false);
        if (_3357_ladder_break)
        {
            break;
        }
        i = _5392 + 1u;
        _5395 = _3527;
        _5392 = i;
        continue;
    }
    RetainedTriple _5394;
    _5394 = _5395;
    RetainedTriple _4275;
    uint i_1;
    uint _5393 = 0u;
    for (;;)
    {
        bool _3373_ladder_break = false;
        do
        {
            if (!(_5393 < 3u))
            {
                _3373_ladder_break = true;
                break;
            }
            bool _4236 = _5393 == 0u;
            uint _5403;
            if (_4236)
            {
                _5403 = a.lo;
            }
            else
            {
                uint _5404;
                if (_5393 == 1u)
                {
                    _5404 = a.hi;
                }
                else
                {
                    _5404 = a.tail;
                }
                _5403 = _5404;
            }
            uint _5406;
            if (_4236)
            {
                _5406 = b.lo;
            }
            else
            {
                uint _5407;
                if (_5393 == 1u)
                {
                    _5407 = b.tail;
                }
                else
                {
                    _5407 = b.hi;
                }
                _5406 = _5407;
            }
            RetainedTriple _3336 = _5394;
            _4275 = RTAccumulateRoundedProduct(_3336, _5403, _5406);
            break;
        } while(false);
        if (_3373_ladder_break)
        {
            break;
        }
        i_1 = _5393 + 1u;
        _5394 = _4275;
        _5393 = i_1;
        continue;
    }
    if (_5394.valid == 0u)
    {
        return _5394;
    }
    uint _4416 = RPAbs(b.tail);
    uint _4420 = RPAbs(a.tail);
    return RTFinishAccumulation(_5394.hi, _5394.lo, _5394.tail, RPUpAdd(RPUpAdd(_5394.error, RPUpAdd(RPUpAdd(RPUpMultiply(RPAbs(a.lo), _4416), RPUpMultiply(_4420, RPAbs(b.lo))), RPUpMultiply(_4420, _4416))), RPUpAdd(RPUpAdd(RPUpMultiply(a.error, RTCenterUpper(b)), RPUpMultiply(b.error, RTCenterUpper(a))), RPUpMultiply(a.error, b.error))));
}

static inline __attribute__((always_inline))
RetainedTriple RTAdd(thread const RetainedTriple& a, thread const RetainedTriple& b)
{
    bool _5461;
    if (a.valid == 0u)
    {
        _5461 = true;
    }
    else
    {
        _5461 = b.valid == 0u;
    }
    if (_5461)
    {
        return RTInvalid();
    }
    uint _4520 = RPZero();
    if (RTExact(a, _4520))
    {
        return b;
    }
    if (RTExact(b, _4520))
    {
        return a;
    }
    RetainedTriple _5463;
    _5463 = a;
    RetainedTriple _4557;
    uint i;
    uint _5462 = 0u;
    for (;;)
    {
        bool _4489_ladder_break = false;
        do
        {
            if (!(_5462 < 3u))
            {
                _4489_ladder_break = true;
                break;
            }
            uint _5467;
            if (_5462 == 0u)
            {
                _5467 = b.hi;
            }
            else
            {
                uint _5468;
                if (_5462 == 1u)
                {
                    _5468 = b.lo;
                }
                else
                {
                    _5468 = b.tail;
                }
                _5467 = _5468;
            }
            RetainedTriple _4477 = _5463;
            _4557 = RTAccumulateTerm(_4477, _5467);
            break;
        } while(false);
        if (_4489_ladder_break)
        {
            break;
        }
        i = _5462 + 1u;
        _5463 = _4557;
        _5462 = i;
        continue;
    }
    RetainedTriple _5464;
    if (_5463.valid == 0u)
    {
        _5464 = _5463;
    }
    else
    {
        _5464 = RTFinishAccumulation(_5463.hi, _5463.lo, _5463.tail, RPUpAdd(_5463.error, b.error));
    }
    return _5464;
}

static inline __attribute__((always_inline))
RetainedTriple RTCenter(thread const RetainedTriple& a)
{
    RetainedTriple _4595 = a;
    _4595.error = RPZero();
    return _4595;
}

static inline __attribute__((always_inline))
RetainedTriple RTEvaluateArithmetic(thread const RetainedTriple& a, thread const RetainedTriple& b, uint op)
{
    bool _5185;
    if (op == 2u)
    {
        _5185 = true;
    }
    else
    {
        _5185 = op == 3u;
    }
    bool direct_product = op == 4u;
    bool square_root = op == 6u;
    RetainedTriple _2077 = RTInvalid();
    uint _2078 = RPZero();
    uint _5206;
    RetainedTriple _5210;
    RetainedTriple _5225;
    if (_5185)
    {
        RetainedTriple _5192;
        if (op == 3u)
        {
            _5192 = RTNegate(b);
        }
        else
        {
            _5192 = b;
        }
        _5225 = _5192;
        _5210 = _2077;
        _5206 = _2078;
    }
    else
    {
        uint _5207;
        RetainedTriple _5211;
        if (direct_product)
        {
            _5211 = b;
            _5207 = _2078;
        }
        else
        {
            uint _5208;
            RetainedTriple _5212;
            if (square_root)
            {
                if (a.valid == 0u)
                {
                    return a;
                }
                bool _5188;
                if (RPScalarEqual(a.hi, _2078))
                {
                    _5188 = RPScalarEqual(a.lo, _2078);
                }
                else
                {
                    _5188 = false;
                }
                bool _5189;
                if (_5188)
                {
                    _5189 = RPScalarEqual(a.tail, _2078);
                }
                else
                {
                    _5189 = false;
                }
                bool _5190;
                if (_5189)
                {
                    _5190 = RPScalarEqual(a.error, _2078);
                }
                else
                {
                    _5190 = false;
                }
                if (_5190)
                {
                    return RTConstant(_2078);
                }
                bool _5191;
                if (RPScalarLessEqual(a.hi, _2078))
                {
                    _5191 = true;
                }
                else
                {
                    _5191 = RPScalarLessEqual(RTAbsoluteLower(a), _2078);
                }
                if (_5191)
                {
                    return _2077;
                }
                _5212 = RTConstant(RPSqrtSeed(a.hi));
                _5208 = _2078;
            }
            else
            {
                bool _5186;
                if (a.valid == 0u)
                {
                    _5186 = true;
                }
                else
                {
                    _5186 = b.valid == 0u;
                }
                if (_5186)
                {
                    return _2077;
                }
                uint _2886 = RTAbsoluteLower(b);
                bool _5187;
                if (RPScalarLessEqual(_2886, _2078))
                {
                    _5187 = true;
                }
                else
                {
                    _5187 = RPScalarGreater(RPRatioUpper(RTAbsoluteUpper(a), _2886), RPFromBits(1577058304u));
                }
                if (_5187)
                {
                    return _2077;
                }
                if (RTExact(b, RPOne()))
                {
                    return a;
                }
                _5212 = RTConstant(RPDivideSeed(a.hi, b.hi));
                _5208 = _2886;
            }
            _5211 = _5212;
            _5207 = _5208;
        }
        _5225 = b;
        _5210 = _5211;
        _5206 = _5207;
    }
    RetainedTriple _5209;
    RetainedTriple _5243;
    _5243 = _2077;
    _5209 = _5210;
    uint phase;
    uint _5277;
    RetainedTriple _5294;
    RetainedTriple _5303;
    uint _5193 = 0u;
    uint _5205 = _5206;
    for (;;)
    {
        bool _1999_ladder_break = false;
        do
        {
            if (!(_5193 < 7u))
            {
                _1999_ladder_break = true;
                break;
            }
            bool correcting = (_5193 & 1u) != 0u;
            RetainedTriple _5244;
            RetainedTriple _5245;
            if (correcting)
            {
                _5277 = _5205;
                _5245 = _5243;
                _5244 = _5209;
            }
            else
            {
                RetainedTriple _5237;
                uint _5238;
                if (_5185)
                {
                    _5238 = _5205;
                    _5237 = _5225;
                }
                else
                {
                    uint _5240;
                    if (square_root)
                    {
                        uint _5241;
                        if (_5193 == 6u)
                        {
                            RetainedTriple _1943 = _5209;
                            uint _3252 = RTAbsoluteLower(_1943);
                            bool _5218;
                            if (_5209.valid == 0u)
                            {
                                _5218 = true;
                            }
                            else
                            {
                                _5218 = RPScalarLessEqual(_5209.hi, _2078);
                            }
                            bool _5219;
                            if (_5218)
                            {
                                _5219 = true;
                            }
                            else
                            {
                                _5219 = RPScalarLessEqual(_3252, _2078);
                            }
                            if (_5219)
                            {
                                return _2077;
                            }
                            _5241 = _3252;
                        }
                        else
                        {
                            bool _5213;
                            if (_5209.valid == 0u)
                            {
                                _5213 = true;
                            }
                            else
                            {
                                _5213 = RPScalarLessEqual(_5209.hi, _2078);
                            }
                            bool _5215;
                            if (_5213)
                            {
                                _5215 = true;
                            }
                            else
                            {
                                _5215 = !RPBounded(_5209.hi, 1560281088u);
                            }
                            if (_5215)
                            {
                                return _2077;
                            }
                            _5241 = _5205;
                        }
                        _5240 = _5241;
                    }
                    else
                    {
                        _5240 = _5205;
                    }
                    RetainedTriple _5233;
                    if (square_root)
                    {
                        _5233 = _5209;
                    }
                    else
                    {
                        RetainedTriple _5234;
                        if (direct_product)
                        {
                            _5234 = a;
                        }
                        else
                        {
                            _5234 = _5225;
                        }
                        _5233 = _5234;
                    }
                    RetainedTriple _1944 = _5233;
                    RetainedTriple _1945 = _5209;
                    RetainedTriple _3324 = RTMultiply(_1944, _1945);
                    if (direct_product)
                    {
                        return _3324;
                    }
                    RetainedTriple _1946 = _3324;
                    _5238 = _5240;
                    _5237 = RTNegate(_1946);
                }
                _5277 = _5238;
                _5245 = _5237;
                _5244 = a;
            }
            RetainedTriple _1947 = _5244;
            RetainedTriple _1948 = _5245;
            RetainedTriple _4468 = RTAdd(_1947, _1948);
            if (_5185)
            {
                return _4468;
            }
            if (correcting)
            {
                RetainedTriple _1949 = _4468;
                RetainedTriple _4591 = RTCenter(_1949);
                bool _5278;
                if (!square_root)
                {
                    _5278 = _4591.valid == 0u;
                }
                else
                {
                    _5278 = false;
                }
                if (_5278)
                {
                    return _4591;
                }
                _5303 = _5243;
                _5294 = _4591;
            }
            else
            {
                if (_5193 == 6u)
                {
                    if (_4468.valid == 0u)
                    {
                        return _4468;
                    }
                    RetainedTriple _1950 = _4468;
                    return RTPack(_5209.hi, _5209.lo, _5209.tail, RPRatioUpper(RTAbsoluteUpper(_1950), _5277));
                }
                uint _5272;
                if (square_root)
                {
                    _5272 = OrderedAdd32(_5209.hi, _5209.hi);
                }
                else
                {
                    _5272 = _5225.hi;
                }
                uint _5273;
                if (!square_root)
                {
                    _5273 = RPAbs(_5272);
                }
                else
                {
                    _5273 = _5272;
                }
                bool _5274;
                if (_4468.valid == 0u)
                {
                    _5274 = true;
                }
                else
                {
                    _5274 = RPScalarGreater(RPRatioUpper(RPAbs(_4468.hi), _5273), RPFromBits(1560281088u));
                }
                if (_5274)
                {
                    return _2077;
                }
                _5303 = RTConstant(RPDivideSeed(_4468.hi, _5272));
                _5294 = _5209;
            }
            break;
        } while(false);
        if (_1999_ladder_break)
        {
            break;
        }
        phase = _5193 + 1u;
        _5243 = _5303;
        _5209 = _5294;
        _5205 = _5277;
        _5193 = phase;
        continue;
    }
    return _2077;
}

static inline __attribute__((always_inline))
RetainedTriple RTNegate_1(thread const RetainedTriple& a)
{
    RetainedTriple _4693 = a;
    _4693.hi = RPFromBits(RPBits(_4693.hi) ^ 2147483648u);
    _4693.lo = RPFromBits(RPBits(_4693.lo) ^ 2147483648u);
    _4693.tail = RPFromBits(RPBits(_4693.tail) ^ 2147483648u);
    return _4693;
}

static inline __attribute__((always_inline))
uint RTCenterUpper_1(thread const RetainedTriple& a)
{
    return RPUpAdd(RPUpAdd(RPAbs(a.hi), RPAbs(a.lo)), RPAbs(a.tail));
}

static inline __attribute__((always_inline))
uint RTAbsoluteUpper_1(thread const RetainedTriple& a)
{
    return RPUpAdd(RTCenterUpper_1(a), a.error);
}

static inline __attribute__((always_inline))
uint RPScalarMax(uint a, uint b)
{
    uint _5471;
    if (PB32Less(a, b))
    {
        _5471 = b;
    }
    else
    {
        _5471 = a;
    }
    return _5471;
}

static inline __attribute__((always_inline))
void ScratchWrite(uint index, uint component, uint rowBase, uint prefix, uint value, threadgroup spvUnsafeArray<uint, 1810>& retainedScratch)
{
    retainedScratch[ScratchAddress(index, component, rowBase, prefix)] = value;
}

static inline __attribute__((always_inline))
void WriteValue(uint index, thread const RetainedTriple& r, uint rowBase, threadgroup spvUnsafeArray<uint, 1810>& retainedScratch)
{
    ScratchWrite(index, 0u, rowBase, 204u, RPBits(r.hi), retainedScratch);
    ScratchWrite(index, 1u, rowBase, 204u, RPBits(r.lo), retainedScratch);
    ScratchWrite(index, 2u, rowBase, 204u, RPBits(r.tail), retainedScratch);
    ScratchWrite(index, 3u, rowBase, 204u, RPBits(r.error), retainedScratch);
    ScratchWrite(index, 4u, rowBase, 204u, r.valid, retainedScratch);
}

static inline __attribute__((always_inline))
bool EvaluateOne(uint index, uint row, uint rowBase, uint programBase, uint registers, device StructuredBuffer& inputs, constant uint& inputsBufferSize, threadgroup spvUnsafeArray<uint, 1810>& retainedScratch)
{
    uint address = (programBase + 43u) + (5u * index);
    uint _1771 = address + 1u;
    uint _1773 = inputs._member0[_1771];
    uint _1774 = address + 2u;
    uint _1777 = address + 3u;
    uint _1780 = address + 4u;
    bool _5163;
    if (_1773 >= registers)
    {
        _5163 = true;
    }
    else
    {
        _5163 = inputs._member0[address] > 11u;
    }
    if (_5163)
    {
        return false;
    }
    RetainedTriple _1794 = RTInvalid();
    RetainedTriple _5173;
    if (inputs._member0[address] == 0u)
    {
        _5173 = RTConstant(RPFromBits(inputs._member0[_1774]));
    }
    else
    {
        RetainedTriple _5174;
        if (inputs._member0[address] == 1u)
        {
            if (inputs._member0[_1774] >= 45u)
            {
                return false;
            }
            _5174 = ReadOriginal(inputs._member0[_1774], row, inputs, inputsBufferSize);
        }
        else
        {
            if (inputs._member0[_1774] >= registers)
            {
                return false;
            }
            RetainedTriple av = ReadValue(inputs._member0[_1774], rowBase, retainedScratch);
            bool _5164;
            if (inputs._member0[address] >= 2u)
            {
                _5164 = inputs._member0[address] <= 6u;
            }
            else
            {
                _5164 = false;
            }
            RetainedTriple _5175;
            if (_5164)
            {
                RetainedTriple _5172;
                if (inputs._member0[address] != 6u)
                {
                    if (inputs._member0[_1777] >= registers)
                    {
                        return false;
                    }
                    _5172 = ReadValue(inputs._member0[_1777], rowBase, retainedScratch);
                }
                else
                {
                    _5172 = _1794;
                }
                RetainedTriple _1684 = av;
                RetainedTriple _1685 = _5172;
                _5175 = RTEvaluateArithmetic(_1684, _1685, inputs._member0[address]);
            }
            else
            {
                RetainedTriple _5176;
                if (inputs._member0[address] == 7u)
                {
                    RetainedTriple _1686 = av;
                    _5176 = RTNegate_1(_1686);
                }
                else
                {
                    RetainedTriple _5177;
                    if (inputs._member0[address] == 8u)
                    {
                        RetainedTriple _1687 = av;
                        uint _4721 = RTAbsoluteUpper_1(_1687);
                        bool _5171;
                        if (av.valid != 0u)
                        {
                            _5171 = RPScalarGreater(_4721, RPZero());
                        }
                        else
                        {
                            _5171 = false;
                        }
                        RetainedTriple _5178;
                        if (_5171)
                        {
                            _5178 = RTConstant(RPPowerOfTwo(RPFloorExponent(_4721) + 1));
                        }
                        else
                        {
                            _5178 = _1794;
                        }
                        _5177 = _5178;
                    }
                    else
                    {
                        if (inputs._member0[_1777] >= registers)
                        {
                            return false;
                        }
                        RetainedTriple bv = ReadValue(inputs._member0[_1777], rowBase, retainedScratch);
                        RetainedTriple _5179;
                        if (inputs._member0[address] == 9u)
                        {
                            bool _5170;
                            if (av.valid != 0u)
                            {
                                _5170 = bv.valid != 0u;
                            }
                            else
                            {
                                _5170 = false;
                            }
                            RetainedTriple _5180;
                            if (_5170)
                            {
                                RetainedTriple _1688 = av;
                                RetainedTriple _1689 = bv;
                                _5180 = RTConstant(RPScalarMax(RTAbsoluteUpper_1(_1688), RTAbsoluteUpper_1(_1689)));
                            }
                            else
                            {
                                _5180 = _1794;
                            }
                            _5179 = _5180;
                        }
                        else
                        {
                            bool _4810 = inputs._member0[address] == 10u;
                            bool _5165;
                            if (_4810)
                            {
                                _5165 = true;
                            }
                            else
                            {
                                _5165 = inputs._member0[address] == 11u;
                            }
                            RetainedTriple _5181;
                            if (_5165)
                            {
                                if (inputs._member0[_1780] >= registers)
                                {
                                    return false;
                                }
                                RetainedTriple cv = ReadValue(inputs._member0[_1780], rowBase, retainedScratch);
                                RetainedTriple _5182;
                                if (cv.valid != 0u)
                                {
                                    bool _5168;
                                    if (_4810)
                                    {
                                        _5168 = RPScalarGreaterEqual(cv.hi, RPZero());
                                    }
                                    else
                                    {
                                        uint _4834 = RPZero();
                                        bool _5166;
                                        if (RPScalarEqual(cv.hi, _4834))
                                        {
                                            _5166 = RPScalarEqual(cv.lo, _4834);
                                        }
                                        else
                                        {
                                            _5166 = false;
                                        }
                                        bool _5167;
                                        if (_5166)
                                        {
                                            _5167 = RPScalarEqual(cv.tail, _4834);
                                        }
                                        else
                                        {
                                            _5167 = false;
                                        }
                                        bool _5169;
                                        if (_5167)
                                        {
                                            _5169 = RPScalarEqual(cv.error, _4834);
                                        }
                                        else
                                        {
                                            _5169 = false;
                                        }
                                        _5168 = _5169;
                                    }
                                    RetainedTriple _5183;
                                    if (_5168)
                                    {
                                        _5183 = av;
                                    }
                                    else
                                    {
                                        _5183 = bv;
                                    }
                                    _5182 = _5183;
                                }
                                else
                                {
                                    _5182 = _1794;
                                }
                                _5181 = _5182;
                            }
                            else
                            {
                                _5181 = _1794;
                            }
                            _5179 = _5181;
                        }
                        _5177 = _5179;
                    }
                    _5176 = _5177;
                }
                _5175 = _5176;
            }
            _5174 = _5175;
        }
        _5173 = _5174;
    }
    RetainedTriple _1690 = _5173;
    WriteValue(_1773, _1690, rowBase, retainedScratch);
    return true;
}

static inline __attribute__((always_inline))
void WriteStored_1(uint base, thread const RetainedTriple& r, device RWStructuredBuffer& outputs, constant uint& outputsBufferSize)
{
    outputs._member0[base] = RPBits(r.hi);
    outputs._member0[base + 1u] = RPBits(r.lo);
    outputs._member0[base + 2u] = RPBits(r.tail);
    outputs._member0[base + 3u] = RPBits(r.error);
    outputs._member0[base + 4u] = r.valid;
}

static inline __attribute__((always_inline))
bool EvaluateProgram(uint row, uint rowBase, uint programBase, uint count, uint registers, uint required, uint lane, device StructuredBuffer& inputs, constant uint& inputsBufferSize, device RWStructuredBuffer& outputs, constant uint& outputsBufferSize, threadgroup uint& retainedStatus, threadgroup spvUnsafeArray<uint, 1810>& retainedScratch)
{
    bool _5144;
    if (count == 0u)
    {
        _5144 = true;
    }
    else
    {
        _5144 = count > count;
    }
    if (_5144)
    {
        return false;
    }
    bool _1597 = lane == 0u;
    if (_1597)
    {
        retainedStatus = 0u;
    }
    uint layerBase = (programBase + 43u) + (5u * count);
    bool _5145;
    if (inputs._member0[layerBase] != 0u)
    {
        _5145 = true;
    }
    else
    {
        _5145 = inputs._member0[layerBase + 73u] != count;
    }
    if (_5145)
    {
        return false;
    }
    uint _1632;
    uint layer;
    uint _5146 = 0u;
    uint _5147 = 0u;
    for (;;)
    {
        bool _1524_ladder_break = false;
        do
        {
            if (!(_5146 < 73u))
            {
                _1524_ladder_break = true;
                break;
            }
            uint _1627 = layerBase + _5146;
            uint _1629 = inputs._member0[_1627];
            _1632 = inputs._member0[_1627 + 1u];
            if (_1629 >= count)
            {
                _1524_ladder_break = true;
                break;
            }
            bool _5148;
            if (_1629 != _5147)
            {
                _5148 = true;
            }
            else
            {
                _5148 = _1632 <= _1629;
            }
            bool _5149;
            if (_5148)
            {
                _5149 = true;
            }
            else
            {
                _5149 = _1632 > count;
            }
            bool _5150;
            if (_5149)
            {
                _5150 = true;
            }
            else
            {
                _5150 = (_1632 - _1629) > 64u;
            }
            if (_5150)
            {
                RejectRetainedProgram(retainedStatus);
            }
            else
            {
                uint node;
                uint _5151 = _1629;
                for (;;)
                {
                    bool _1542_ladder_break = false;
                    do
                    {
                        if (!(_5151 < _1632))
                        {
                            _1542_ladder_break = true;
                            break;
                        }
                        bool _1673 = EvaluateOne(_5151, row, rowBase, programBase, registers, inputs, inputsBufferSize, retainedScratch);
                        if (!_1673)
                        {
                            RejectRetainedProgram(retainedStatus);
                        }
                        break;
                    } while(false);
                    if (_1542_ladder_break)
                    {
                        break;
                    }
                    node = _5151 + 1u;
                    _5151 = node;
                    continue;
                }
            }
            if (retainedStatus != 0u)
            {
                return false;
            }
            break;
        } while(false);
        if (_1524_ladder_break)
        {
            break;
        }
        layer = _5146 + 1u;
        _5147 = _1632;
        _5146 = layer;
        continue;
    }
    uint i;
    uint _5159 = 0u;
    for (;;)
    {
        bool _1562_ladder_break = false;
        do
        {
            if (!(_5159 < required))
            {
                _1562_ladder_break = true;
                break;
            }
            uint _4950 = (programBase + 3u) + _5159;
            bool _5161;
            if (inputs._member0[_4950] >= registers)
            {
                _5161 = true;
            }
            else
            {
                _5161 = ReadValue(inputs._member0[_4950], rowBase, retainedScratch).valid == 0u;
            }
            if (_5161)
            {
                return false;
            }
            break;
        } while(false);
        if (_1562_ladder_break)
        {
            break;
        }
        i = _5159 + 1u;
        _5159 = i;
        continue;
    }
    if (_1597)
    {
        uint i_1;
        uint _5160 = 0u;
        for (;;)
        {
            bool _1578_ladder_break = false;
            do
            {
                if (!(_5160 < 40u))
                {
                    _1578_ladder_break = true;
                    break;
                }
                RetainedTriple _1509 = ReadValue(inputs._member0[(programBase + 3u) + _5160], rowBase, retainedScratch);
                WriteStored_1((rowBase + 4u) + (_5160 * 5u), _1509, outputs, outputsBufferSize);
                break;
            } while(false);
            if (_1578_ladder_break)
            {
                break;
            }
            i_1 = _5160 + 1u;
            _5160 = i_1;
            continue;
        }
    }
    return true;
}

kernel void main0(constant spvDescriptorSetBuffer0& spvDescriptorSet0 [[buffer(0)]], constant uint* spvBufferSizeConstants [[buffer(29)]], uint3 gl_WorkGroupID [[threadgroup_position_in_grid]], uint gl_LocalInvocationIndex [[thread_index_in_threadgroup]])
{
    threadgroup uint retainedStatus;
    threadgroup spvUnsafeArray<uint, 1810> retainedScratch;
    constant uint& spvDescriptorSet0_inputsBufferSize = spvDescriptorSet0.spvBufferSizeConstants[1];
    constant uint& spvDescriptorSet0_outputsBufferSize = spvDescriptorSet0.spvBufferSizeConstants[2];
    uint2 _131 = uint2((spvDescriptorSet0_inputsBufferSize - 0) / 4, 4u);
    uint ni = _131.x;
    uint2 _137 = uint2((spvDescriptorSet0_outputsBufferSize - 0) / 4, 4u);
    uint no = _137.x;
    bool _5041;
    if (ni < 1u)
    {
        _5041 = true;
    }
    else
    {
        _5041 = gl_WorkGroupID.y != 0u;
    }
    bool _5042;
    if (_5041)
    {
        _5042 = true;
    }
    else
    {
        _5042 = gl_WorkGroupID.z != 0u;
    }
    if (_5042)
    {
        return;
    }
    uint _170 = (*spvDescriptorSet0.inputs)._member0[0];
    bool _5043;
    if (_170 == 0u)
    {
        _5043 = true;
    }
    else
    {
        _5043 = _170 > 65536u;
    }
    bool _5044;
    if (_5043)
    {
        _5044 = true;
    }
    else
    {
        _5044 = gl_WorkGroupID.x >= _170;
    }
    bool _5045;
    if (_5044)
    {
        _5045 = true;
    }
    else
    {
        _5045 = ni < ((1u + (225u * _170)) + 43u);
    }
    bool _5046;
    if (_5045)
    {
        _5046 = true;
    }
    else
    {
        _5046 = (no % _170) != 0u;
    }
    if (_5046)
    {
        return;
    }
    uint programBase = 1u + (225u * _170);
    uint _214 = (*spvDescriptorSet0.inputs)._member0[programBase];
    uint _215 = programBase + 1u;
    uint _217 = (*spvDescriptorSet0.inputs)._member0[_215];
    uint rowWords = no / _170;
    uint rowBase = gl_WorkGroupID.x * rowWords;
    if (rowWords < 204u)
    {
        return;
    }
    uint i;
    uint _5047 = gl_LocalInvocationIndex;
    for (;;)
    {
        bool _49_ladder_break = false;
        do
        {
            if (!(_5047 < rowWords))
            {
                _49_ladder_break = true;
                break;
            }
            (*spvDescriptorSet0.outputs)._member0[rowBase + _5047] = 0u;
            break;
        } while(false);
        if (_49_ladder_break)
        {
            break;
        }
        i = _5047 + 1u;
        _5047 = i;
        continue;
    }
    bool _5048;
    if ((*spvDescriptorSet0.inputs)._member0[programBase + 2u] != 40u)
    {
        _5048 = true;
    }
    else
    {
        _5048 = _214 > 65536u;
    }
    bool _5049;
    if (_5048)
    {
        _5049 = true;
    }
    else
    {
        _5049 = _217 != 362u;
    }
    bool _5050;
    if (_5049)
    {
        _5050 = true;
    }
    else
    {
        _5050 = ni != ((((programBase + 43u) + (5u * _214)) + 73u) + 1u);
    }
    bool _5051;
    if (_5050)
    {
        _5051 = true;
    }
    else
    {
        _5051 = rowWords != (204u + (5u * _217));
    }
    if (_5051)
    {
        return;
    }
    uint i_1;
    uint _5052 = 0u;
    for (;;)
    {
        bool _73_ladder_break = false;
        do
        {
            if (!(_5052 < 45u))
            {
                _73_ladder_break = true;
                break;
            }
            if (ReadOriginal(_5052, gl_WorkGroupID.x, (*spvDescriptorSet0.inputs), spvDescriptorSet0_inputsBufferSize).valid == 0u)
            {
                return;
            }
            break;
        } while(false);
        if (_73_ladder_break)
        {
            break;
        }
        i_1 = _5052 + 1u;
        _5052 = i_1;
        continue;
    }
    RetainedTriple chart = ReadOriginal(44u, gl_WorkGroupID.x, (*spvDescriptorSet0.inputs), spvDescriptorSet0_inputsBufferSize);
    uint _1133 = RPOne();
    RetainedTriple _13 = chart;
    bool _5053;
    if (!RTExact(_13, _1133))
    {
        RetainedTriple _14 = chart;
        _5053 = !RTExact(_14, RPScalarNegate(_1133));
    }
    else
    {
        _5053 = false;
    }
    if (_5053)
    {
        return;
    }
    RetainedTriple _15 = ReadOriginal(1u, gl_WorkGroupID.x, (*spvDescriptorSet0.inputs), spvDescriptorSet0_inputsBufferSize);
    uint _1345 = RPZero();
    bool _5054;
    if (!RPScalarEqual(RTAbsoluteUpper(_15), _1345))
    {
        RetainedTriple _16 = ReadOriginal(3u, gl_WorkGroupID.x, (*spvDescriptorSet0.inputs), spvDescriptorSet0_inputsBufferSize);
        _5054 = !RPScalarEqual(RTAbsoluteUpper(_16), _1345);
    }
    else
    {
        _5054 = false;
    }
    if (_5054)
    {
        return;
    }
    RetainedTriple _17 = ReadOriginal(0u, gl_WorkGroupID.x, (*spvDescriptorSet0.inputs), spvDescriptorSet0_inputsBufferSize);
    bool _5055;
    if (RTExact(_17, _1345))
    {
        RetainedTriple _18 = ReadOriginal(2u, gl_WorkGroupID.x, (*spvDescriptorSet0.inputs), spvDescriptorSet0_inputsBufferSize);
        _5055 = RTExact(_18, _1345);
    }
    else
    {
        _5055 = false;
    }
    bool _5056;
    if (_5055)
    {
        RetainedTriple _19 = ReadOriginal(3u, gl_WorkGroupID.x, (*spvDescriptorSet0.inputs), spvDescriptorSet0_inputsBufferSize);
        _5056 = RTExact(_19, _1345);
    }
    else
    {
        _5056 = false;
    }
    if (_5056)
    {
        if (gl_LocalInvocationIndex != 0u)
        {
            return;
        }
        uint i_2;
        uint _5057 = 0u;
        for (;;)
        {
            bool _104_ladder_break = false;
            do
            {
                if (!(_5057 < 40u))
                {
                    _104_ladder_break = true;
                    break;
                }
                RetainedTriple value = ReadOriginal(4u + _5057, gl_WorkGroupID.x, (*spvDescriptorSet0.inputs), spvDescriptorSet0_inputsBufferSize);
                if (gl_LocalInvocationIndex == 0u)
                {
                    RetainedTriple _5058;
                    if ((_5057 % 8u) == 4u)
                    {
                        RetainedTriple _21 = value;
                        _5058 = RTNegate(_21);
                    }
                    else
                    {
                        _5058 = value;
                    }
                    RetainedTriple _22 = _5058;
                    WriteStored((rowBase + 4u) + (_5057 * 5u), _22, (*spvDescriptorSet0.outputs), spvDescriptorSet0_outputsBufferSize);
                }
                break;
            } while(false);
            if (_104_ladder_break)
            {
                break;
            }
            i_2 = _5057 + 1u;
            _5057 = i_2;
            continue;
        }
        if (gl_LocalInvocationIndex == 0u)
        {
            (*spvDescriptorSet0.outputs)._member0[rowBase] = 1u;
        }
        return;
    }
    bool _1492 = EvaluateProgram(gl_WorkGroupID.x, rowBase, programBase, _214, _217, 40u, gl_LocalInvocationIndex, (*spvDescriptorSet0.inputs), spvDescriptorSet0_inputsBufferSize, (*spvDescriptorSet0.outputs), spvDescriptorSet0_outputsBufferSize, retainedStatus, retainedScratch);
    if (!_1492)
    {
        return;
    }
    if (gl_LocalInvocationIndex == 0u)
    {
        (*spvDescriptorSet0.outputs)._member0[rowBase] = 1u;
    }
}

