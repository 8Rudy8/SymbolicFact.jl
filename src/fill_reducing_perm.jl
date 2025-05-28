"""
    fill_reducing_perm(A)

Return another matrix with permuted rows to reduce the fill-in
during factorization.
"""
function fill_reducing_perm(A)
	if(issymmetric(A))
		perm, iperm = Metis.permutation(A)
		B = A[perm,perm]
	else
		perm, iperm = Metis.permutation(A'*A)
		B = A[:,perm]
	end
	return B
end