#include "sirius/kernels/portable_binary32.h"
#include "word_io.h"
#include <iostream>
int main(int argc, char** argv) {
    try {
        if (argc != 4) throw std::runtime_error("usage: cpu_probe input expected actual");
        const auto input=ReadWords(argv[1]),expected=ReadWords(argv[2]);
        CheckShape(input,expected);
        using namespace sirius::portable_binary32;
        std::vector<std::uint32_t> actual(expected.size());
        for (std::size_t row=0;row<input[0];++row) {
            const auto op=input[1+4*row],a=input[2+4*row],b=input[3+4*row];
            auto& high=actual[4*row]; auto& low=actual[4*row+1];auto& valid=actual[4*row+2];valid=1;
            switch(op) {
                case 0:high=PB32Add(a,b);break;
                case 1:high=PB32Subtract(a,b);break;
                case 2:high=PB32Multiply(a,b);break;
                case 3:high=PB32Divide(a,b);break;
                case 4:high=PB32Sqrt(a);break;
                case 5:high=PB32Equal(a,b);break;
                case 6:high=PB32Less(a,b);break;
                case 7:high=PB32LessEqual(a,b);break;
                case 8:high=PB32FromUnsigned(a);break;
                case 9:high=PB32FromSignedBits(a);break;
                case 10:high=PB32Abs(a);break;
                case 11:high=PB32Negate(a);break;
                case 12:high=PB32Classify(a);break;
                case 13:{const auto p=PB32ProductResidual(a,b);high=p.high;low=p.low;valid=p.valid;break;}
                default:valid=0;break;
            }
        }
        std::size_t mismatches=0;
        for(std::size_t i=0;i<actual.size();++i) {
            if(actual[i]==expected[i])continue;
            if(mismatches<8)std::cerr<<"row="<<i/4<<" column="<<i%4<<" expected="<<expected[i]<<" actual="<<actual[i]<<'\n';
            ++mismatches;
        }
        WriteWords(argv[3],actual);
        std::cout<<"{\"records\":"<<input[0]<<",\"words\":"<<actual.size()<<",\"mismatches\":"<<mismatches<<"}\n";
        return mismatches==0?0:1;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
