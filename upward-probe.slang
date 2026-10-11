import retained_scalar;
import retained_pair;
[[vk::binding(0,0)]] StructuredBuffer<uint> inputs;
[[vk::binding(1,0)]] RWStructuredBuffer<uint> outputs;
[shader("compute")] [numthreads(64,1,1)]
void main(uint3 id:SV_DispatchThreadID) {
 uint n,stride,m,other;inputs.GetDimensions(n,stride);outputs.GetDimensions(m,other);
 if(n<1u||id.y!=0u||id.z!=0u)return;
 uint count=inputs[0];if(count>(n-1u)/3u||count>m/5u||id.x>=count)return;
 uint at=1u+3u*id.x,a=inputs[at+1],b=inputs[at+2];
 uint domain=a>=0x00800000u&&a<=0x5e800000u&&b>=0x00800000u&&b<=0x5e800000u&&(a>>23)+(b>>23)>=128u?1u:0u;
 uint out_at=5u*id.x;outputs[out_at]=RPBits(RPUpMultiply(RPFromBits(a),RPFromBits(b)));
 outputs[out_at+1]=domain;outputs[out_at+2]=a;outputs[out_at+3]=b;outputs[out_at+4]=inputs[at];
}
