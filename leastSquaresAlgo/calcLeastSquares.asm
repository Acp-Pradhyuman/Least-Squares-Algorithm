; .386
.model flat, c
.const
public LsEpsilon
LsEpsilon real8 1.0e-12
.code
; calcLeastSquaresASM(const double* x, const double* y, int n, double *m, 
; double *b);
calcLeastSquaresASM proc
    push ebp
    mov ebp, esp
    sub esp, 8      ; creating a local variable (8bytes for double)

    ; stack contents
    ; local variable    -> esp
    ; ebp               -> ebp
    ; ret addr
    ; x                 -> 4bytes
    ; y
    ; n                 -> 4bytes
    ; m
    ; b

    xor eax, eax
    mov ecx, [ebp+16]   ; ecx <- n
    test ecx, ecx
    jle done

    mov eax, [ebp+8]    ; eax <- x
    mov edx, [ebp+12]   ; edx <- y

    fldz                ; sumXX = 0
    fldz                ; sumXY = 0
    fldz                ; sumY  = 0             
    fldz                ; sumX  = 0

    ; contents of x87 stack
    ; sumX
    ; sumY
    ; sumXY
    ; sumXX

@@: fld real8 ptr [eax]
    fld st(0)
    fld st(0)
    fld real8 ptr [edx]

    ; contents of x87 stack
    ; y
    ; x
    ; x
    ; x
    ; sumX
    ; sumY
    ; sumXY
    ; sumXX

    ; note: make sure that the x87 stack contains only upto 8 values
    ; beyond that the stack overflows and we get unexpected output
    ; st(0) to st(7)

    ; Add ST(i) to ST(0) and store result in ST(i)
    fadd st(5), st(0)   ; st(5) = sumY + y

    ; Multiply ST(1) by ST(0), store result in ST(1), 
    ; and pop the register stack.
    fmulp               ;  y*x

    ; contents of x87 stack
    ; y*x
    ; x
    ; x
    ; sumX
    ; sumY + y
    ; sumXY
    ; sumXX

    ; Add ST(0) to ST(i), store result in ST(i), and pop the register stack
    faddp st(5), st(0)  ; st(6) = sumXY + y*x, and pop

    ; contents of x87 stack
    ; x
    ; x
    ; sumX
    ; sumY + y
    ; sumXY + y*x
    ; sumXX

    fadd st(2), st(0)  ; st(2) = sumX + x, and pop

    ; contents of x87 stack
    ; x
    ; x
    ; sumX + x
    ; sumY + y
    ; sumXY + y*x
    ; sumXX

    fmulp   ; st(1) = st(1) * st(0), and pop

    ; contents of x87 stack
    ; x*x
    ; sumX + x
    ; sumY + y
    ; sumXY + y*x
    ; sumXX

    faddp st(4), st(0)  ; sumXX + x*x

    add eax, 8
    add edx, 8
    dec ecx
    jnz @B

; double denom = n * sumXX - sumX * sumX;
    ; FILD m32int -> Push m32int onto the FPU register stack.
    fild dword ptr [ebp+16]     ; st(0) <- n

    ; contents of x87 stack
    ; n
    ; sumX (final value)
    ; sumY (final value)
    ; sumXY (final value)
    ; sumXX (final value)
    
    ; FMUL ST(i), ST(0) -> Multiply ST(i) by ST(0) & store result in ST(i)
    ; FMUL ST(0), ST(i) -> Multiply ST(0) by ST(i) and store result in ST(0).
    fmul st(0), st(4)   ; n * sumXX

    ; contents of x87 stack
    ; n*sumXX
    ; sumX
    ; sumY
    ; sumXY
    ; sumXX

    fld st(1)
    fld st(0)

    ; contents of x87 stack
    ; sumX
    ; sumX
    ; n*sumXX
    ; sumX
    ; sumY
    ; sumXY
    ; sumXX

    fmulp

    ; Subtract ST(0) from ST(1), store result in ST(1), & pop register stack.
    fsubp

    ; FST m64fp -> Copy ST(0) to m64fp.
    fst real8 ptr [ebp-8]   ; local variable (stores denom)

    ; contents of x87 stack
    ; n*sumXX - sumX*sumX = denom
    ; sumX
    ; sumY
    ; sumXY
    ; sumXX

    ; fabs is a floating-point instruction that converts the value in ST(0) 
    ; to its absolute value
    fabs
    fld real8 ptr [LsEpsilon]

    ; contents of x87 stack
    ; LsEpsilon
    ; abs(n*sumXX - sumX*sumX) = abs(denom)
    ; sumX
    ; sumY
    ; sumXY
    ; sumXX

    ; Compare ST(0) with ST(i), set status flags accordingly, 
    ; and pop register stack.
    fcomip st(0), st(1)

    ; contents of x87 stack
    ; abs(n*sumXX - sumX*sumX) = abs(denom)
    ; sumX
    ; sumY
    ; sumXY
    ; sumXX

    ; FSTP ST(i) -> Copy ST(0) to ST(i) and pop register stack.
    fstp st(0)

    ; LsEpsilon >= fabs(denom)
    jae invalidDenom    ; jump if epsilon >= fabs(denom)

    ; contents of x87 stack
    ; sumX
    ; sumY
    ; sumXY
    ; sumXX
    

    ; *m = (n * sumXY - sumX * sumY) / denom; (slope)
    fild dword ptr [ebp+16]     ; st(0) <- n
    fmul st(0), st(3)           ; n*sumXY
    fld st(2)                   ; sumY
    fld st(2)                   ; sumX
    fmulp                       ; sumX*sumY
    fsubp                       ; n*sumXY - sumX*sumY

    ; contents of x87 stack
    ; n*sumXY - sumX*sumY
    ; sumX
    ; sumY
    ; sumXY
    ; sumXX

    ; FDIV m64fp -> Divide ST(0) by m64fp and store result in ST(0).
    ; ebp-8 -> local variable (contains denom)
    fdiv real8 ptr [ebp-8]  
    mov eax, [ebp+20]

    ; contents of x87 stack
    ; (n*sumXY - sumX*sumY) / denom
    ; sumX
    ; sumY
    ; sumXY
    ; sumXX   

    fstp real8 ptr [eax]        ; *m = (n * sumXY - sumX * sumY) / denom

    ; contents of x87 stack
    ; sumX
    ; sumY
    ; sumXY
    ; sumXX 

    ; *b = (sumXX*sumY - sumX*sumXY) / denom; (intercept)
    ; FXCH ST(i) -> Exchange the contents of ST(0) and ST(i).
    fxch st(3)

    ; contents of x87 stack
    ; sumXX
    ; sumY
    ; sumXY
    ; sumX

    fmulp
    fxch st(2)
    fmulp

    ; contents of x87 stack
    ; sumX*sumXY
    ; sumXX*sumY

    fsubp

    ; contents of x87 stack
    ; sumXX*sumY - sumX*sumXY

    fdiv real8 ptr [ebp-8]      ; local variable (denom value)

    ; *b = (sumXX*sumY - sumX*sumXY) / denom;

    ; contents of x87 stack
    ; (sumXX*sumY - sumX*sumXY) / denom

    mov eax, [ebp+24]
    fstp real8 ptr [eax] 
    mov eax, 1

done:
    mov esp, ebp
    pop ebp
    ret

invalidDenom:
    fstp st(0)
    fstp st(0)
    fstp st(0)
    fstp st(0)

    xor eax, eax
    mov esp, ebp
    pop ebp
    ret

calcLeastSquaresASM endp

; avxAddPackedFp(double x[], double y[], int size, double* sumX, double* sumY);
avxAddPackedFp proc
    push ebp
    mov ebp, esp

    ; contents of stack
    ; ebp       -> ebp (4-bytes)
    ; ret addr  (4-bytes)
    ; x         (4-bytes)
    ; y         (4-bytes)
    ; size      (4-bytes)
    ; sumX      (4-bytes)
    ; sumY

    mov eax, [ebp+8]    ; eax <- x
    mov edx, [ebp+12]   ; edx <- y
    mov ecx, [ebp+16]   ; ecx <- size
    shr ecx, 3
    ; dec ecx

    vzeroall

@@:
    ; VMOVAPD ymm1, ymm2/m256
    ; Move aligned packed double precision floating-point values from 
    ; ymm2/mem to ymm1
    vmovapd ymm2, ymmword ptr [eax]
    vaddpd ymm1, ymm1, ymm2
    vmovapd ymm0, ymmword ptr [eax+32]
    vaddpd ymm1, ymm1, ymm0

    vmovapd ymm5, ymmword ptr [edx]
    vaddpd ymm4, ymm4, ymm5
    vmovapd ymm3, ymmword ptr [edx+32]
    vaddpd ymm4, ymm4, ymm3
    add eax, 64
    add edx, 64
    loop @B

    ; VHADDPD (VEX.256 Encoded Version)
    ; DEST[63:0] := SRC1[127:64] + SRC1[63:0]
    ; DEST[127:64] := SRC2[127:64] + SRC2[63:0]
    ; DEST[191:128] := SRC1[255:192] + SRC1[191:128]
    ; DEST[255:192] := SRC2[255:192] + SRC2[191:128]
    ; if ymm1 have [y3, y2, y1, y0]
    ; where y3, y2, y1, y0 are double data results
    ; After vhaddpd, 
    ; YMM1 now contains as qwords y3+y2 (two times), 
    ; y1+y0 (two times)
    ; ymm1 have [y3+y2, y3+y2, y1+y0, y1+y0]
    vhaddpd	ymm1, ymm1, ymm1
    vhaddpd	ymm4, ymm4, ymm4
    ; ymm2 = [y3+y2, y3+y2, y1+y0, y1+y0]
    vmovapd ymm2, ymm1
    vmovapd ymm5, ymm4

    ; VPERMPD ymm1, ymm2/m256, imm8
    ; DEST[63:0] := (SRC[255:0] >> (IMM8[1:0] * 64))[63:0];
    ; DEST[127:64] := (SRC[255:0] >> (IMM8[3:2] * 64))[63:0];
    ; DEST[191:128] := (SRC[255:0] >> (IMM8[5:4] * 64))[63:0];
    ; DEST[255:192] := (SRC[255:0] >> (IMM8[7:6] * 64))[63:0];
    ; DEST[MAXVL-1:256] := 0
    ; 2 = 00000010
    ; IMM8[1:0] = 10 → Shift by 2 * 64 = 128 bits.
    ; IMM8[3:2] = 00 → Shift by 0 * 64 = 0 bits.
    ; IMM8[5:4] = 00 → Shift by 0 * 64 = 0 bits.
    ; IMM8[7:6] = 00 → Shift by 0 * 64 = 0 bits.
    ; so ymm2 = [y3+y2, y3+y2, y1+y0, y1+y0]
    ; becomes ymm2 = [y1+y0, y1+y0, y1+y0, y3+y2]
    vpermpd ymm2, ymm2, 2
    ; vpermpd ymm2, ymm2, 3 also works
    vpermpd ymm5, ymm5, 2

    ; ymm2 = [y1+y0, y1+y0, y1+y0, y3+y2]
    ; ymm1 = [y3+y2, y3+y2, y1+y0, y1+y0]
    ; ymm0 = [x    , x    , x    , y3+y2+y1+y0]
    vaddpd ymm0, ymm1, ymm2
    vaddpd ymm3, ymm4, ymm5

    ; vmovapd real8 ptr [edx], ymm0
    ; vmovsd qword ptr [edx], xmm0
    mov eax, [ebp+20]   ; eax <- sumX 
    mov edx, [ebp+24]   ; edx <- sumY
    vmovsd qword ptr [eax], xmm0
    vmovsd qword ptr [edx], xmm3

    vzeroupper

    pop ebp
    ret
avxAddPackedFp endp

; fmaPackedFp(double x[], double y[], int size, double* sumXX, double* sumXY);
fmaPackedFp proc
    push ebp
    mov ebp, esp

    ; contents of stack
    ; ebp       -> ebp (4-bytes)
    ; ret addr  (4-bytes)
    ; x         (4-bytes)
    ; y         (4-bytes)
    ; size      (4-bytes)
    ; sumXX     (4-bytes)
    ; sumXY

    mov eax, [ebp+8]    ; eax <- x
    mov edi, [ebp+12]   ; edi <- y
    mov ecx, [ebp+16]   ; ecx <- size
    shr ecx, 2

    xor esi, esi

    ; remember to clear ymm registers since it might have previous results
    vzeroall

@@:
    ; VMOVAPD ymm1, ymm2/m256
    ; Move aligned packed double precision floating-point values from 
    ; ymm2/mem to ymm1
    vmovapd ymm1, ymmword ptr [eax + esi]   ; ymm1 <- x
    vmovapd ymm2, ymmword ptr [edi + esi]   ; ymm2 <- y

    ; VFMADD231PD ymm1, ymm2, ymm3/m256
    vfmadd231pd ymm3, ymm1, ymm1            ; ymm3 calculates sumXX
    vfmadd231pd ymm4, ymm1, ymm2            ; ymm4 calculates sumXY
    add esi, 32
    loop @B

    vhaddpd	ymm3, ymm3, ymm3
    vmovapd ymm2, ymm3
    vpermpd ymm2, ymm2, 2
    vaddpd ymm0, ymm2, ymm3

    vhaddpd	ymm4, ymm4, ymm4
    vmovapd ymm2, ymm4
    vpermpd ymm2, ymm2, 2
    vaddpd ymm1, ymm2, ymm4

    mov eax, [ebp+20]   ; eax <- sumX 
    mov edx, [ebp+24]   ; edx <- sumY

    ; vmovapd real8 ptr [edx], ymm0
    vmovsd qword ptr [eax], xmm0
    vmovsd qword ptr [edx], xmm1

    vzeroupper

    pop ebp
    ret
fmaPackedFp endp
end