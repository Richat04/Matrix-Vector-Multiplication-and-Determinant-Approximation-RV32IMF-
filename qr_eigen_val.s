# the data section starts from here. all the inputs are fed here
.data
# N - just a placeholder for the value6 which is the length of the matrix
N:      .word 6
# input matrix A = 6X6 matrix
A:      .float -109.0, 3.0, 3.0, 4.0, 5.0, 6.0        
        .float 2.0, 4.0, 2.0, 1.2, 2.9, 1.0
        .float 3.7, 0.0, 1.0, 2.4, 21.0, 0.0
        .float 1.0, 2.0, 3.5, 4.0, 3.0, 2.2
        .float 2.0, -9.0, 20.0, 0.9, 2.0, 1.0
        .float 1.8, 3.0, 2.0, 4.0, 1.0, 5.0
#input matrix X - 6X1 matrix
X:      .float 11.0, 2.52, -3.87, 4.0, -5.96, 0.0
#output Y - 6X1 matrix. 24 bits space is given since 6X4 bytes
Y:      .space 24

det_result: .float 0.0
det_qr_result: .float 0.0
comparison_flag: .word 0
# N,A,X are our inputs while Y,det_result, det_qr_result and comparison flag are our outputs
A_temp: .space 144   
# A_temp is the copied matrix used for guassian elimination
# Working matrices for QR method
work_matrix: .space 144
Q: .space 144
R: .space 144
temp_matrix: .space 144
temp_col: .space 24

# Constants for QR
CONST_ONE: .float 1.0
EPS: .float 0.001
MAX_ITER: .word 100
# these constants are used for comparison checks
# Output messages
#msq_titel	.asciz "///////////////Results///////////"
msg_y:          .asciz "Vector Y = A*X:\n"
msg_space:      .asciz " "
msg_newline:    .asciz "\n"
msg_sum:        .asciz "Sum of Y: "
msg_det_gauss:  .asciz "Determinant (Gaussian Elimination): "
msg_det_qr:     .asciz "Determinant (QR Eigenvalue Method): "
msg_comp_flag:  .asciz "Comparison Flag (sum(Y) < det): "
#msq_mat_A	.asciz " Matrix A is....."
# main function ( the actual computation and printing) starts from here
.text
.globl main
#.main
main:
    # Initialize constants
    jal ra, init_constants  #it will jump to the next function which is init_constants

    # matvec: Y = A * X
    #loading all the input in the registers
    la a0, A
    la a1, X
    la a2, Y
    lw a3, N
    jal ra, matvec_mul  #jusp to the function matvec_mul for the matrix multiplicatio i.e A*X = Y

    # determinant (Gaussian elimination)
    la a0, A
    lw a1, N
    jal ra, det_comp   #jump to det_comp - computation of gaussian determinant starts from this function

    # store determinant to memory
    la t0, det_result          # load det_result in t0
    fsw fa0, 0(t0)  #store the computed determinnt (t0) to floating point reg fa0

    # determinant (QR eigenvalue method)
    la a0, A
    lw a1, N
    jal ra, det_qr  #jumps to det_qr for computing determinant using qr method

    # store QR determinant
    la t0, det_qr_result
    fsw fa0, 0(t0)   #same as above, first load the reuslt i t0 and then store it in fa0

    # compare sum(Y) with det (Gaussian)
    la a0, Y
    lw a1, N
    la t0, det_result
    flw fa1, 0(t0)
    jal ra, compare_sum_det  #jump to comare_sum_det function

    la t0, comparison_flag #load comparison flag in t0
    sw a0, 0(t0)  #store the result in a0

    # Print results
    jal ra, print_results # finally jump to print results

    # exit
    li a7, 10
    ecall  #function ends
#///////////////////////////#
init_constants:
    la t0, CONST_ONE  #just load some constant values in float registrors and return
    flw f31, 0(t0)      # f31 = 1.0
    ret
#////////// matrix multiplication///////////#
matvec_mul:  #matrix multiplication starts from here
    mv t2, a3          # t2 = N    
    li t0, 0           # row index t0(i) = 0

matvec_outer:
    bge t0, t2, matvec_done #if row index is equal to 6 then jump to matvec_done function

    # ft0 = 0.0  (accumulator)
    fcvt.s.w ft0, zero   #convert int 0 to float (single precision
	#it is an accumulator i.e. at the end it has one element of matrix Y
    li t1, 0           # j = 0 - column index

matvec_inner:
    bge t1, t2, matvec_store #if column index(j) == 6, move to matvec_store (store the result)

    # addr = A + (i*N + j) * 4
    mul t3, t0, t2
    add t3, t3, t1
    slli t3, t3, 2
    add t3, a0, t3
    flw ft1, 0(t3)
    # the above 5 lines are written to compute the linear index, in short find the location of A[i][j]
	# then load the value of A[i][j] in ft1
    # X[j]
    slli t3, t1, 2
    add t3, a1, t3
    flw ft2, 0(t3)
	#did the same for X matrix
    # ft0 += ft1 * ft2
    fmul.s ft1, ft1, ft2
    fadd.s ft0, ft0, ft1
	#multiply A[i][j] with X[j] and add it to the accumulator
    addi t1, t1, 1  #increase j by one
    j matvec_inner  #loop continues

matvec_store:
    slli t3, t0, 2
    add t3, a2, t3
    fsw ft0, 0(t3)
	#compute the address for Y i.e. Y + i*4 are store the accumulated values there
    addi t0, t0, 1  #increase row count
    j matvec_outer #loop continues

matvec_done:
    ret  #when i = 6, return

#/////////////Gaussian//////////
det_comp:
    
    addi sp, sp, -40  #reserve 40 bytes in stack pointer
    sw s0, 0(sp)
    sw s1, 4(sp)
    sw s2, 8(sp)
    sw s3, 12(sp)
    sw s4, 16(sp)
    sw s5, 20(sp)
    sw s6, 24(sp)
    sw s7, 28(sp)
    sw ra, 32(sp)
	#saved all the registers
    mv s0, a0      # s0 = A base
    mv s1, a1      # s1 = N
    la s2, A_temp  # s2 = A_temp base

    # copy A -> A_temp (N*N floats)
    mul t0, s1, s1   # t0 = N*N
    li t1, 0
copy_loop:
    bge t1, t0, copy_done
    slli t2, t1, 2
    add t3, s0, t2   
    #add t2, t2, t0
    #slli t2, t2, 2
    #add t2, s2, t2   
    add t4, s2, t2
    flw ft0, 0(t3)
    fsw ft0, 0(t4)
    addi t1, t1, 1
    j copy_loop
    #copy everything in A_temp so that we dont change original A
copy_done:

    # set determinant = 1.0 in ft5
    li t0, 0x3f800000  # IEEE 754 representation of 1.0
    fmv.w.x ft5, t0 #this reg will be multiplied by diagonal ele and negated on row swaps

    li t0, 0      # k = 0 (pivot column)

elim_outer:
    bge t0, s1, elim_after   # finished elimination (when k = 6)

    # load pivot A_temp[k][k] into ft0
    mul t1, t0, s1
    add t1, t1, t0
    slli t1, t1, 2
    add t1, s2, t1
    flw ft0, 0(t1)

    # Check if pivot == 0.0
    fcvt.s.w ft6, zero    # ft6 = 0.0
    feq.s t3, ft0, ft6    # t3 = 1 if pivot == 0.0

    beq t3, zero, pivot_nonzero   # if pivot != 0 continue

    # pivot == 0.0 -> try to find a row r > k with non-zero in column k
    addi s3, t0, 1     # s3 = r (start r = k+1)
find_nonzero: #now serachig for thr row r>k with non zero in column k
    bge s3, s1, pivot_all_zero  # no non-zero pivot found
    # load A_temp[r][k]
    mul t5, s3, s1
    add t5, t5, t0
    slli t5, t5, 2
    add t5, s2, t5
    flw ft1, 0(t5)
    feq.s t6, ft1, ft6   # t6 = 1 if A_temp[r][k] == 0.0
    beq t6, zero, found_row_swap #t6 = 0 means that ft1 is non zero
    addi s3, s3, 1 #increase r
    j find_nonzero

found_row_swap:
    # swap rows k and r (for j = 0..N-1)
    li s4, 0      # s4 = j = 0
row_swap_loop:
    bge s4, s1, row_swap_done
    # addr_k = A_temp[k][j]
    mul t4, t0, s1
    add t4, t4, s4
    slli t4, t4, 2
    add t4, s2, t4
    # addr_r = A_temp[r][j]
    mul t5, s3, s1
    add t5, t5, s4
    slli t5, t5, 2
    add t5, s2, t5

    flw ft1, 0(t4)
    flw ft2, 0(t5)
    fsw ft2, 0(t4)
    fsw ft1, 0(t5)
	#take out all elements of r one by one and swap it with those of row k
    addi s4, s4, 1
    j row_swap_loop
row_swap_done:
    # flip determinant sign: ft5 = ft5 * -1.0
    li t4, 0xbf800000  # IEEE 754 for -1.0
    fmv.w.x ft7, t4
    fmul.s ft5, ft5, ft7
	#sign flipped due to row swap
    # reload pivot ft0 (now from swapped row)
    mul t1, t0, s1
    add t1, t1, t0
    slli t1, t1, 2
    add t1, s2, t1
    flw ft0, 0(t1)

    j pivot_continue

pivot_all_zero:
    # pivot column all zeros -> determinant = 0
    fcvt.s.w fa0, zero
    # restore saved regs
    lw s0, 0(sp)
    lw s1, 4(sp)
    lw s2, 8(sp)
    lw s3, 12(sp)
    lw s4, 16(sp)
    lw s5, 20(sp)
    lw s6, 24(sp)
    lw s7, 28(sp)
    lw ra, 32(sp)
    addi sp, sp, 40
    ret
	#return early
pivot_nonzero:  #if pivot is non_zero then continue further
pivot_continue:
    # For rows i = k+1 .. N-1
    addi t1, t0, 1

elim_middle:
    bge t1, s1, elim_next_k

    # multiplier = A_temp[i][k] / pivot (ft2)
    mul t2, t1, s1
    add t2, t2, t0
    slli t2, t2, 2
    add t2, s2, t2
    flw ft1, 0(t2)      # ft1 = A_temp[i][k]
    fdiv.s ft2, ft1, ft0

    # inner loop over j = k .. N-1
    mv t3, t0

elim_inner:
    bge t3, s1, elim_next_i

    # A_temp[k][j] -> ft3
    mul t4, t0, s1
    add t4, t4, t3
    slli t4, t4, 2
    add t4, s2, t4
    flw ft3, 0(t4)

    # A_temp[i][j] -> ft4
    mul t5, t1, s1
    add t5, t5, t3
    slli t5, t5, 2
    add t5, s2, t5
    flw ft4, 0(t5)

    # ft4 = ft4 - multiplier * ft3
    fmul.s ft6, ft2, ft3
    fsub.s ft4, ft4, ft6
    fsw ft4, 0(t5)
	#basically performing row subtraction
    addi t3, t3, 1
    #la a0, msg_newline
    #li a7, 4
    #ecall
    j elim_inner  

elim_next_i:
    addi t1, t1, 1
    j elim_middle  #move to next i

elim_next_k:
    addi t0, t0, 1
    j elim_outer  #after finishing this for k, move to next pivot

elim_after:
    # multiply diagonal elements to get determinant
    li t0, 0

det_mult_loop:
    bge t0, s1, det_mult_done

    mul t1, t0, s1
    add t1, t1, t0
    slli t1, t1, 2
    add t1, s2, t1
    flw ft0, 0(t1)
	#loop over index 0 - n-1. Multiply all the diagonal elements to get determinant
    fmul.s ft5, ft5, ft0
    addi t0, t0, 1
    j det_mult_loop

det_mult_done:
    fmv.s fa0, ft5
	#move value to fa0
    # restore saved regs & return
    lw s0, 0(sp)
    lw s1, 4(sp)
    lw s2, 8(sp)
    lw s3, 12(sp)
    lw s4, 16(sp)
    lw s5, 20(sp)
    lw s6, 24(sp)
    lw s7, 28(sp)
    lw ra, 32(sp)
    addi sp, sp, 40
    #la a0, msg_newline
    #li a7, 4
    #ecall
    ret

#/////////////QR METHOD/////////
det_qr:
    addi sp, sp, -16
    sw s0, 12(sp)
    sw s1, 8(sp)
    sw ra, 4(sp)   #allocate 16 bytes on stack and save registers s0,s1,ra

    mv s0, a0
    mv s1, a1  #these are local reg for this part, move values of A and N to s0 and s1

    # Copy A to work_matrix
    la a1, work_matrix
    mv a2, s1
    jal ra, copy_matrix  #jump to copy_matrix
    la a0, work_matrix
    mv a1, s1
    jal ra, qr_iteration  # jump to qr_iteration
    la a0, work_matrix
    mv a1, s1
    jal ra, diag_product #jump to diag_product to compute determinant by product of diagonals
    lw s0, 12(sp)
    lw s1, 8(sp)
    lw ra, 4(sp) #restore registers
    addi sp, sp, 16
    ret #return

copy_matrix:
    mul t0, a2, a2 # a2 = N
    li t1, 0 #t0 contains n*n and t1 copies each value
copy_mat_loop:
    bge t1, t0, copy_mat_done # when t1 = n*n, the loop terminates
    slli t2, t1, 2
    add t3, a0, t2
    flw ft0, 0(t3)
    add t4, a1, t2
    fsw ft0, 0(t4)
    addi t1, t1, 1  #increment t1 by 1 and the loop is executed again
    j copy_mat_loop
copy_mat_done: #to terminante copy loop program
    ret

qr_iteration:
    addi sp, sp, -32
    sw ra, 28(sp)
    sw s0, 24(sp)
    sw s1, 20(sp)
    sw s2, 16(sp)  #save registers
    
    mv s0, a0
    mv s1, a1
    
    lw s2, MAX_ITER #load max_iter
    li t0, 0  #initialise t0 to 0
    
qr_iter_loop:
    bge t0, s2, qr_iter_done  #if max_ietration reached then terminate
    
    addi sp, sp, -4 #save iteration on stack
    sw t0, 0(sp)
    
    mv a0, s0
    mv a1, s1
    jal ra, qr_decompose #jump to qr decompose
    
    la a0, R # load r
    la a1, Q # load Q
    la a2, temp_matrix
    mv a3, s1
    jal ra, matrix_multiply #multiply R and Q
    
    la a0, temp_matrix
    mv a1, s0
    mv a2, s1
    jal ra, copy_matrix #copy temp_matrxi to work_matrix
    
    mv a0, s0
    mv a1, s1
    jal ra, check_convergence #compares off diagonal matrxi to EPS and returns 1 if converged in a0
    
    lw t0, 0(sp)
    addi sp, sp, 4
    
    bnez a0, qr_iter_done #if converged, jump to qr_iter_done
    
    addi t0, t0, 1 # else increment t0
    j qr_iter_loop #loop continues
    
qr_iter_done:
    lw ra, 28(sp)
    lw s0, 24(sp)
    lw s1, 20(sp)
    lw s2, 16(sp)
    addi sp, sp, 32 #now the work_matrix is approx upper traingualr with eigen values
    ret


qr_decompose:
    addi sp, sp, -32  #now the gram-schmidt QR starts
    sw ra, 28(sp)
    sw s0, 24(sp)
    sw s1, 20(sp)  # same as above
    sw s2, 16(sp)
    
    mv s0, a0
    mv s1, a1 #same as above
    
    la a0, Q
    mv a1, s1
    jal ra, make_identity  #initialise matrix Q to identity
    
    la a0, R
    mv a1, s1
    jal ra, zero_matrix  #this zeroes R
    
    li s2, 0
qr_col_loop:
    bge s2, s1, qr_decomp_done #s2 is current column index (k from 0 - n-1)
    
    li t0, 0
copy_col_loop:
    bge t0, s1, gram_schmidt #jump to gram-schmidt if cloumn index>= total rows
    mul t1, t0, s1
    add t1, t1, s2 #add column offset
    slli t1, t1, 2
    add t1, s0, t1 #address of matrix element
    flw ft0, 0(t1) #load the value
    
    la t2, temp_col  #load address of temp_col
    slli t3, t0, 2
    add t3, t2, t3
    fsw ft0, 0(t3)
    
    addi t0, t0, 1  #increase t0 by 1
    j copy_col_loop
    
gram_schmidt:
    li t0, 0
orth_loop:
    bge t0, s2, compute_norm_qr # jump to norm
    
    li t1, 0
    fcvt.s.w ft0, t1  #inititalise dot product to 0
    li t1, 0
    
dot_qr_loop:  #compute dot product between Q[0:t0] and temp_col
    bge t1, s1, subtract_proj    
    la t2, Q
    mul t3, t1, s1
    add t3, t3, t0
    slli t3, t3, 2    
    #add t2, t2, t0
    #slli t2, t2, 2
    #add t2, s2, t2    
    add t3, t2, t3
    flw ft1, 0(t3) #load q element    
    la t2, temp_col
    slli t3, t1, 2
    add t3, t2, t3
    flw ft2, 0(t3) #load temp_col elemnt    
    fmul.s ft3, ft1, ft2
    fadd.s ft0, ft0, ft3 #compute dot product and sum it
    
    addi t1, t1, 1 #increment t1
    j dot_qr_loop # jump to this
    
subtract_proj:
    la t1, R #load address of R
    mul t2, t0, s1
    add t2, t2, s2
    slli t2, t2, 2
    add t2, t1, t2
    fsw ft0, 0(t2) #store the dot product
    
    li t1, 0
sub_proj_loop:
    bge t1, s1, next_orth
    
    la t2, temp_col
    slli t3, t1, 2
    add t3, t2, t3
    flw ft1, 0(t3) #this is temp_col[t1]
    
    la t2, Q #load q in t2
    
    #mul t4, t2, s3
    #add t4, t4, t1
    
    mul t4, t1, s1
    add t4, t4, t0
    slli t4, t4, 2
    add t4, t2, t4
    flw ft2, 0(t4) #load from memory   
    fmul.s ft3, ft0, ft2
    fsub.s ft1, ft1, ft3 # temp_col[t1] -= R[t0,s2] * Q[t1, t0]   
    la t2, temp_col #load temp_col
    slli t3, t1, 2
    add t3, t2, t3
    fsw ft1, 0(t3)
    
    addi t1, t1, 1
    j sub_proj_loop
    
next_orth:
    addi t0, t0, 1 # repeat until columns 0 - s2-1 are processed
    j orth_loop
    
compute_norm_qr:
    li t0, 0
    fcvt.s.w ft0, t0 #convert integer to float
    li t0, 0
    
norm_qr_loop:
    bge t0, s1, normalize_qr
    
    la t1, temp_col
    slli t2, t0, 2
    add t2, t1, t2
    flw ft1, 0(t2)
    fmul.s ft2, ft1, ft1 #compute square
    fadd.s ft0, ft0, ft2 #sum of squaers
    
    addi t0, t0, 1 #increment t0
    j norm_qr_loop
    
normalize_qr:
    fsqrt.s ft0, ft0 #compute square root
    
    la t0, R #load the address of R
    mul t1, s2, s1
    add t1, t1, s2
    slli t1, t1, 2
    add t1, t0, t1
    fsw ft0, 0(t1) #store this element
    
    li t0, 0
store_q_col:
    bge t0, s1, next_qr_col  #if all rows done then jump to this function
    
    la t1, temp_col
    slli t2, t0, 2 #sfiht left
    add t2, t1, t2
    flw ft1, 0(t2) #load from mem
    fdiv.s ft1, ft1, ft0 
    
    la t1, Q #load address of Q
    mul t2, t0, s1
    add t2, t2, s2
    slli t2, t2, 2
    add t2, t1, t2
    fsw ft1, 0(t2) #store in menory
    
    addi t0, t0, 1 #increase index by 1
    j store_q_col #repest
    
next_qr_col:
    addi s2, s2, 1 #move to next colum
    j qr_col_loop #repeat
    
qr_decomp_done:
    lw ra, 28(sp)
    lw s0, 24(sp)
    lw s1, 20(sp)
    lw s2, 16(sp) #restore the saved reg
    addi sp, sp, 32
    ret #fucntion done, retuen

matrix_multiply:
    addi sp, sp, -32  #create soace in stack
    sw s0, 28(sp)
    sw s1, 24(sp)
    sw s2, 20(sp)
    sw s3, 16(sp)   
    mv s0, a0 #move addresses to these reg
    mv s1, a1
    mv s2, a2
    mv s3, a3    
    li t0, 0 #row_index = 0
mul_i_loop:
    bge t0, s3, mul_done #done, complete   
    li t1, 0
mul_j_loop:
    bge t1, s3, mul_next_i    # go to next row if this satisfies
    li t2, 0
    fcvt.s.w ft0, t2    #convert ot float, accumulator is 0 
    li t2, 0
mul_k_loop:
    bge t2, s3, mul_store #store result if k>=6
    
    mul t3, t0, s3
    add t3, t3, t2
    slli t3, t3, 2
    add t3, s0, t3
    flw ft1, 0(t3) #load value of matrix A
    
    mul t4, t2, s3
    add t4, t4, t1
    #add t4, s1, t4
    slli t4, t4, 2
    add t4, s1, t4
    flw ft2, 0(t4) #load value from amtrix B
    
    fmul.s ft3, ft1, ft2 #multiply
    fadd.s ft0, ft0, ft3 #add to accumulate
    
    addi t2, t2, 1 #increment
    j mul_k_loop #again
    
mul_store:
    mul t3, t0, s3
    
    #mul t4, t2, s3
    #add t4, t4, t1
    
    add t3, t3, t1
    slli t3, t3, 2
    add t3, s2, t3
    fsw ft0, 0(t3) #store in C[i][j]    
    addi t1, t1, 1 # increment
    j mul_j_loop #jump
    
mul_next_i:
    addi t0, t0, 1 #increment
    j mul_i_loop #jump
    
mul_done:
    lw s0, 28(sp)
    lw s1, 24(sp)
    lw s2, 20(sp) #load words
    lw s3, 16(sp)
    addi sp, sp, 32
    ret #finally khatam


check_convergence:
    mv t0, a0 #address of matrxi
    mv t1, a1
    
    la t2, EPS #small epsilon for threshhold
    flw ft0, 0(t2) #load value
    
    li t2, 1
conv_loop:
    bge t2, t1, converged #if all rows are done then it mena s converged  
    addi t3, t2, -1
    mul t4, t2, t1
    add t4, t4, t3
    slli t4, t4, 2 #shift left by 2
    add t4, t0, t4
    flw ft1, 0(t4) #load subdiagonal elemts
    
    fabs.s ft1, ft1
    flt.s t4, ft0, ft1 #compare with espisoln
    bnez t4, not_converged #if greater then not converged
    
    addi t2, t2, 1
    j conv_loop
    
converged:
    li a0, 1 #ret 1 if converegd
    ret
    
not_converged:
    li a0, 0 #if not convereged then reutrn 0
    ret

diag_product: #computing the diagonal product = dtereminant mil gaya
    addi sp, sp, -16
    sw s0, 12(sp)
    sw s1, 8(sp)
    sw s2, 4(sp)   
    mv s0, a0 # address of matrix
    mv s1, a1    #load 6
    fmv.s fa0, f31    
    li s2, 0 #index to 0
diag_loop:
    bge s2, s1, diag_done #if greater than or equal to 6 then all diagonals done
    
    addi t0, s2, 1
    bge t0, s1, single_element
    
    mul t1, t0, s1
    add t1, t1, s2
    slli t1, t1, 2 #shift left
    add t1, s0, t1
    flw ft0, 0(t1) #load values from A matrix = A[i+1][i]
    
    fabs.s ft0, ft0 #set to the absolute value
    la t1, EPS #load epsilon
    flw ft1, 0(t1)
    flt.s t1, ft0, ft1
    bnez t1, single_element #branch to this if near zero   
    mul t1, s2, s1
    add t1, t1, s2
    slli t1, t1, 2 #shift
    add t1, s0, t1
    flw ft0, 0(t1) #load A[i][i]   
    addi t2, s2, 1
    mul t1, t2, s1    
    #la a0, msg_newline
    #li a7, 4  
    add t1, t1, t2
    slli t1, t1, 2
    add t1, s0, t1
    flw ft1, 0(t1) #load A[i+1][i+1]
    
    fmul.s ft2, ft0, ft1    
    mul t1, s2, s1
    addi t2, s2, 1
    add t1, t1, t2
    slli t1, t1, 2 #shift left
    add t1, s0, t1
    flw ft0, 0(t1) #load from mem   
    addi t2, s2, 1
    mul t1, t2, s1
    #la a0, msg_newline
    #li a7, 4  
    add t1, t1, s2
    slli t1, t1, 2 #shift left
    add t1, s0, t1
    flw ft1, 0(t1) #load from memory    
    fmul.s ft3, ft0, ft1    
    fsub.s ft2, ft2, ft3    
    fmul.s fa0, fa0, ft2    
    addi s2, s2, 2
    j diag_loop #jump to this label
    
single_element:
    mul t3, s2, s1
    add t3, t3, s2
    slli t3, t3, 2 #shift left
    add t3, s0, t3
    flw ft0, 0(t3) #load from emoery   
    fmul.s fa0, fa0, ft0   
    addi s2, s2, 1
    j diag_loop #jump to thid label
    
diag_done:
    lw s0, 12(sp)
    lw s1, 8(sp)
    lw s2, 4(sp) #loading words
    addi sp, sp, 16 
    ret # return

make_identity:
    li t0, 0
    fcvt.s.w ft0, t0 #convert int to float
    
    li t0, 0
identity_i:
    bge t0, a1, identity_done #if t1 is greater or equal then jump to this fuction
    li t1, 0
identity_j:
    bge t1, a1, identity_next_i #if t1 is greater or equal then jump to this fucnotion
    
    mul t2, t0, a1
    add t2, t2, t1
    slli t2, t2, 2
    add t2, a0, t2
    
    bne t0, t1, identity_zero
    fsw f31, 0(t2) #store to memoery
    j identity_next_j #again jump
identity_zero:
    fsw ft0, 0(t2) #store in memoery
identity_next_j:
    addi t1, t1, 1 #again increment
    j identity_j #again jump
identity_next_i:
    addi t0, t0, 1 #increment t0
    j identity_i #jump to identity_i
identity_done:
    ret
# the next lines are for zeroing out n*n floats at memory a0
zero_matrix:
    mul t0, a1, a1
    li t1, 0
    fcvt.s.w ft0, t1 #convert into float
    li t1, 0
zero_loop:
    bge t1, t0, zero_done # when t1 = 6, loop terminates
    slli t2, t1, 2
    add t2, a0, t2
    fsw ft0, 0(t2)
    addi t1, t1, 1 #increment t1
    j zero_loop
zero_done:
    ret #return
#-------compare program-----
compare_sum_det:
    mv t1, a1
    li t0, 0
    fcvt.s.w ft0, zero   # sum = 0.0

compare_loop:
    bge t0, t1, compare_check
    slli t2, t0, 2
    add t2, a0, t2
    flw ft1, 0(t2)
    fadd.s ft0, ft0, ft1
    addi t0, t0, 1
    j compare_loop

compare_check:
    # flt.s rd, fs1, fs2 is rd = 1 if fs1 < fs2
    flt.s t0, ft0, fa1
    mv a0, t0
    ret  #return

#_________print program_________
print_results:
    addi sp, sp, -4
    sw ra, 0(sp)
    # Print "Vector Y = A*X:"
    la a0, msg_y
    li a7, 4
    ecall
    # Print Y vector elements
    la t0, Y
    lw t1, N
    li t2, 0
print_y_loop:
    bge t2, t1, print_y_done
    flw fa0, 0(t0)
    li a7, 2         # print float
    ecall   
    la a0, msg_space
    li a7, 4
    ecall    
    addi t0, t0, 4
    addi t2, t2, 1
    j print_y_loop #jump to print _ loop
print_y_done:
    la a0, msg_newline
    li a7, 4
    ecall

    # Compute and print sum of Y
    la a0, msg_sum
    li a7, 4
    ecall

    la t0, Y #load address of Y
    lw t1, N  #load value of N
    li t2, 0
    fcvt.s.w ft0, zero
print_sum_loop:
    bge t2, t1, print_sum_done #if t2 = 6 move on to printing it
    flw ft1, 0(t0)
    fadd.s ft0, ft0, ft1
    addi t0, t0, 4
    addi t2, t2, 1 #increment t2
    j print_sum_loop
print_sum_done:
    fmv.s fa0, ft0  #move it to fa0 and print sum of Y
    li a7, 2
    ecall
    la a0, msg_newline #print message line
    li a7, 4
    ecall

    # Print determinant (Gaussian)
    #printing label string and the determinant stored in fa0
    la a0, msg_det_gauss
    li a7, 4
    ecall
    la t0, det_result #load det value
    flw fa0, 0(t0)
    li a7, 2
    ecall
    la a0, msg_newline
    li a7, 4
    ecall

    # Print determinant (QR)
    #printing label string and the determinant stored in fa0
    la a0, msg_det_qr
    li a7, 4
    ecall
    la t0, det_qr_result #load det value
    flw fa0, 0(t0)
    li a7, 2
    ecall
    la a0, msg_newline
    li a7, 4  #load message value to print as well as the result
    ecall

    # Print comparison flag
    la a0, msg_comp_flag  #load message values
    li a7, 4
    ecall
    la t0, comparison_flag
    lw a0, 0(t0)
    li a7, 1         # print integer
    ecall
    la a0, msg_newline
    li a7, 4
    ecall

    lw ra, 0(sp)  #restore ra
    addi sp, sp, 4 #adjust sp
    ret