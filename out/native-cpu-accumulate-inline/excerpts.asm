
RTAccumulateTerm_0___un_3C_s_5B_vyRetainedTriple_0_5D__3E_vyf
   ac7d0:	mov    %rdi,%rax
   ac7d3:	vpcmpeqd %ymm2,%ymm2,%ymm2
   ac7d7:	vtestps %ymm2,%ymm1
   ac7dc:	vmovdqu 0x80(%rsi),%ymm3
   ac7e4:	jb     ac828 <RTAccumulateTerm_0___un_3C_s_5B_vyRetainedTriple_0_5D__3E_vyf+0x58>
   ac7e6:	vpxor  %xmm2,%xmm2,%xmm2
   ac7ea:	vpcmpeqd %ymm2,%ymm3,%ymm3
   ad19f:	vpxor  %xmm3,%xmm3,%xmm3
   ad1a3:	vxorps %xmm5,%xmm5,%xmm5
   ad1a7:	vpxor  %xmm6,%xmm6,%xmm6
   ad1ab:	vmovaps %ymm1,0x80(%rax)
   ad1b3:	vmovaps %ymm6,0x60(%rax)
   ad1b8:	vmovaps %ymm5,0x40(%rax)
   ad1bd:	vmovaps %ymm3,0x20(%rax)
   ad1c2:	vmovaps %ymm0,(%rax)
   ad1c6:	vzeroupper
   ad1c9:	ret

RTAccumulateProduct_0___un_3C_s_5B_vyRetainedTriple_0_5D__3E_vyfvyf
   ad1d0:	push   %rbp
   ad1d1:	mov    %rsp,%rbp
   ad1d4:	push   %rbx
   ad1d5:	and    $0xffffffffffffffe0,%rsp
   ad1d9:	sub    $0x500,%rsp
   ad1e0:	vmovdqa %ymm2,%ymm9
   ad1e4:	mov    %rdi,%rbx
   ad49c:	vsubps %ymm3,%ymm2,%ymm2
   ad4a0:	vsubps %ymm2,%ymm1,%ymm1
   ad4a4:	vmovaps %ymm1,0x100(%rsp)
   ad4ad:	lea    0x300(%rsp),%rdi
   ad4b5:	vmovdqa %ymm7,%ymm1
   ad4b9:	vmovaps %ymm9,0xc0(%rsp)
   ad4c2:	call   ac7d0 <RTAccumulateTerm_0___un_3C_s_5B_vyRetainedTriple_0_5D__3E_vyf>
   ad4c7:	vmovaps 0x380(%rsp),%ymm0
   ad4d0:	vmovaps 0x300(%rsp),%ymm1
   ad4d9:	vmovaps 0x320(%rsp),%ymm2
   ad4e2:	vmovaps 0x340(%rsp),%ymm3
   ad4eb:	vmovaps 0x360(%rsp),%ymm4
   ad4f4:	vmovaps %ymm1,0x3a0(%rsp)
   ad4fd:	vmovaps %ymm2,0x3c0(%rsp)
   ad506:	vmovaps %ymm3,0x3e0(%rsp)
   ad50f:	vmovaps %ymm4,0x400(%rsp)
   ad518:	vmovaps %ymm0,0x420(%rsp)
   ad521:	lea    0x260(%rsp),%rdi
   ad529:	lea    0x3a0(%rsp),%rsi
   ad531:	vmovaps 0x100(%rsp),%ymm0
   ad53a:	vmovaps 0xe0(%rsp),%ymm1
   ad543:	call   ac7d0 <RTAccumulateTerm_0___un_3C_s_5B_vyRetainedTriple_0_5D__3E_vyf>
   ad548:	vmovdqa 0x80(%rsp),%ymm13
   ad551:	vmovaps 0x260(%rsp),%ymm0
   ad55a:	vmovaps 0x280(%rsp),%ymm1
   ad563:	vmovaps 0x2a0(%rsp),%ymm2
   ad56c:	vmovdqa 0x2c0(%rsp),%ymm4
   ad575:	vmovdqa 0x2e0(%rsp),%ymm3
   ad57e:	vxorps %xmm5,%xmm5,%xmm5
   ad582:	vpcmpeqd %ymm5,%ymm3,%ymm5
   ad586:	vpandn %ymm13,%ymm5,%ymm5
   ad58b:	vpand  0xc0(%rsp),%ymm5,%ymm6
   ad594:	vmovmskps %ymm6,%eax
   ad598:	test   %eax,%eax
   ad59a:	je     ad7a2 <RTAccumulateProduct_0___un_3C_s_5B_vyRetainedTriple_0_5D__3E_vyfvyf+0x5d2>
   ad5a0:	vpcmpeqd %ymm7,%ymm7,%ymm7
   ad5a4:	vtestps %ymm7,%ymm6
   ad5a9:	vpbroadcastd 0xbcc8f2(%rip),%ymm8        # c79ea4 <_ZN17sirius_cpu_bridge17transport_programE+0xf9a4>
   adbcd:	vmovdqa %ymm7,0x20(%rbx)
   adbd2:	vmovdqa %ymm6,(%rbx)
   adbd6:	mov    %rbx,%rax
   adbd9:	lea    -0x8(%rbp),%rsp
   adbdd:	pop    %rbx
   adbde:	pop    %rbp
   adbdf:	vzeroupper
   adbe2:	ret
   adbe3:	vbroadcastss 0xbcc290(%rip),%ymm7        # c79e7c <_ZN17sirius_cpu_bridge17transport_programE+0xf97c>
   adbec:	jmp    adb11 <RTAccumulateProduct_0___un_3C_s_5B_vyRetainedTriple_0_5D__3E_vyfvyf+0x941>

RTAccumulateRoundedProduct_0___un_3C_s_5B_vyRetainedTriple_0_5D__3E_vyfvyf
   ade90:	push   %rbp
   ade91:	mov    %rsp,%rbp
   ade94:	push   %rbx
   ade95:	and    $0xffffffffffffffe0,%rsp
   ade99:	sub    $0x280,%rsp
   adea0:	vmovdqa %ymm2,%ymm8
   adea4:	mov    %rdi,%rbx
   ae8aa:	vpcmpeqd %ymm4,%ymm4,%ymm4
   ae8ae:	vpsubd %ymm4,%ymm3,%ymm3
   ae8b2:	vblendvps %ymm1,%ymm6,%ymm3,%ymm6
   ae8b8:	vmovaps %ymm6,(%rsp)
   ae8bd:	lea    0x1c0(%rsp),%rdi
   ae8c5:	vmovaps %ymm2,%ymm1
   ae8c9:	call   ac7d0 <RTAccumulateTerm_0___un_3C_s_5B_vyRetainedTriple_0_5D__3E_vyf>
   ae8ce:	vmovdqa 0x1c0(%rsp),%ymm0
   ae8d7:	vmovaps 0x1e0(%rsp),%ymm1
   ae8e0:	vmovaps 0x200(%rsp),%ymm2
   ae8e9:	vmovdqa 0x220(%rsp),%ymm4
   ae8f2:	vmovdqa 0x240(%rsp),%ymm3
   ae8fb:	vxorps %xmm5,%xmm5,%xmm5
   ae8ff:	vpcmpeqd %ymm5,%ymm3,%ymm5
   ae903:	vmovdqa 0x20(%rsp),%ymm15
   ae909:	vpandn %ymm15,%ymm5,%ymm5
   ae90e:	vpand  0x40(%rsp),%ymm5,%ymm7
   ae914:	vmovmskps %ymm7,%eax
   ae918:	test   %eax,%eax
   ae91a:	je     aedf9 <RTAccumulateRoundedProduct_0___un_3C_s_5B_vyRetainedTriple_0_5D__3E_vyfvyf+0xf69>
   ae920:	vpcmpeqd %ymm8,%ymm8,%ymm8
   ae925:	vtestps %ymm8,%ymm7
   ae92a:	vmovdqa 0x100(%rsp),%ymm10
   aeca3:	vpsubd %ymm3,%ymm2,%ymm2
   aeca7:	vblendvps %ymm1,%ymm5,%ymm2,%ymm5
   aecad:	vmovaps %ymm5,0x20(%rsp)
   aecb3:	lea    0x120(%rsp),%rdi
   aecbb:	vmovaps %ymm15,%ymm1
   aecbf:	vmovaps %ymm15,(%rsp)
   aecc4:	call   ac7d0 <RTAccumulateTerm_0___un_3C_s_5B_vyRetainedTriple_0_5D__3E_vyf>
   aecc9:	vmovdqa (%rsp),%ymm15
   aecce:	vmovdqa 0x120(%rsp),%ymm0
   aecd7:	vmovaps 0x140(%rsp),%ymm1
   aece0:	vmovaps 0x160(%rsp),%ymm2
   aece9:	vmovdqa 0x180(%rsp),%ymm4
   aecf2:	vmovdqa 0x1a0(%rsp),%ymm3
   aecfb:	vxorps %xmm5,%xmm5,%xmm5
   aecff:	vpcmpeqd %ymm5,%ymm3,%ymm5
   aed03:	vpandn %ymm15,%ymm5,%ymm5
   aed08:	vmovmskps %ymm5,%eax
   aed0c:	test   %eax,%eax
   aed0e:	je     aedf9 <RTAccumulateRoundedProduct_0___un_3C_s_5B_vyRetainedTriple_0_5D__3E_vyfvyf+0xf69>
   aed14:	vpcmpeqd %ymm7,%ymm7,%ymm7
   aed18:	vtestps %ymm7,%ymm5
   aed1d:	vmovdqa 0x40(%rsp),%ymm10
   aed23:	vpand  %ymm4,%ymm10,%ymm6
   aef3c:	vmovaps %ymm14,0x60(%rbx)
   aef41:	vmovaps %ymm13,0x40(%rbx)
   aef46:	vmovaps %ymm12,0x20(%rbx)
   aef4b:	vmovaps %ymm11,(%rbx)
   aef4f:	mov    %rbx,%rax
   aef52:	lea    -0x8(%rbp),%rsp
   aef56:	pop    %rbx
   aef57:	pop    %rbp
   aef58:	vzeroupper
   aef5b:	ret

RTAdd_0___un_3C_s_5B_vyRetainedTriple_0_5D__3E_un_3C_s_5B_vyRetainedTriple_0_5D__3E_
   b71e0:	push   %rbp
   b71e1:	mov    %rsp,%rbp
   b71e4:	push   %r14
   b71e6:	push   %rbx
   b71e7:	and    $0xffffffffffffffe0,%rsp
   b71eb:	sub    $0x900,%rsp
   b71f2:	mov    %rdx,%r14
   b7d4a:	vmovaps %ymm10,0xa0(%rsp)
   b7d53:	vmovaps %ymm9,0x20(%rsp)
   b7d59:	vmovaps %ymm0,0xc0(%rsp)
   b7d62:	vmovaps %ymm1,0xe0(%rsp)
   b7d6b:	vmovups (%r14),%ymm0
   b7d70:	lea    0x5e0(%rsp),%rdi
   b7d78:	call   ac7d0 <RTAccumulateTerm_0___un_3C_s_5B_vyRetainedTriple_0_5D__3E_vyf>
   b7d7d:	vmovaps 0x660(%rsp),%ymm0
   b7d86:	vmovaps 0x5e0(%rsp),%ymm1
   b7d8f:	vmovaps 0x600(%rsp),%ymm2
   b7d98:	vmovaps 0x620(%rsp),%ymm3
   b7da1:	vmovaps 0x640(%rsp),%ymm4
   b7daa:	vmovaps %ymm1,0x720(%rsp)
   b7db3:	vmovaps %ymm2,0x740(%rsp)
   b7dbc:	vmovaps %ymm3,0x760(%rsp)
   b7dc5:	vmovaps %ymm4,0x780(%rsp)
   b7dce:	vmovaps %ymm0,0x7a0(%rsp)
   b7dd7:	vmovups 0x20(%r14),%ymm0
   b7ddd:	lea    0x540(%rsp),%rdi
   b7de5:	lea    0x720(%rsp),%rsi
   b7ded:	vmovaps 0xe0(%rsp),%ymm1
   b7df6:	call   ac7d0 <RTAccumulateTerm_0___un_3C_s_5B_vyRetainedTriple_0_5D__3E_vyf>
   b7dfb:	vmovaps 0x5c0(%rsp),%ymm0
   b7e04:	vmovaps 0x540(%rsp),%ymm1
   b7e0d:	vmovaps 0x560(%rsp),%ymm2
   b7e16:	vmovaps 0x580(%rsp),%ymm3
   b7e1f:	vmovaps 0x5a0(%rsp),%ymm4
   b7e28:	vmovaps %ymm1,0x680(%rsp)
   b7e31:	vmovaps %ymm2,0x6a0(%rsp)
   b7e3a:	vmovaps %ymm3,0x6c0(%rsp)
   b7e43:	vmovaps %ymm4,0x6e0(%rsp)
   b7e4c:	vmovaps %ymm0,0x700(%rsp)
   b7e55:	vmovdqu 0x40(%r14),%ymm0
   b7e5b:	lea    0x4a0(%rsp),%rdi
   b7e63:	lea    0x680(%rsp),%rsi
   b7e6b:	vmovaps 0xe0(%rsp),%ymm1
   b7e74:	call   ac7d0 <RTAccumulateTerm_0___un_3C_s_5B_vyRetainedTriple_0_5D__3E_vyf>
   b7e79:	vmovdqa 0x4a0(%rsp),%ymm2
   b8ce9:	vmovaps %ymm12,0x40(%rbx)
   b8cee:	vmovaps %ymm11,0x20(%rbx)
   b8cf3:	vmovaps %ymm10,(%rbx)
   b8cf7:	mov    %rbx,%rax
   b8cfa:	lea    -0x10(%rbp),%rsp
   b8cfe:	pop    %rbx
   b8cff:	pop    %r14
   b8d01:	pop    %rbp
   b8d02:	vzeroupper
   b8d05:	ret
