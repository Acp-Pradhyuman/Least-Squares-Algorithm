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
end