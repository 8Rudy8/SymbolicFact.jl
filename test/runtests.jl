using Test,SparseArrays, MatrixDepot, SymbolicFact
#symmetric matrix
    D = sparse([2,1,2,5,2,4,5,6,7,1,3,7,1,2,3,4,5,6,7,8,9,3,4,4,6,7,7,8,8,8,9,9,9],
    [3,4,4,6,7,7,8,8,8,9,9,9,1,2,3,4,5,6,7,8,9,2,1,2,5,2,4,5,6,7,1,3,7],ones(Int,33))
    #non-symmetricmatrix
    C = sparse([1,2,2,2,3,3,4,4,5,5,6,6,7,7,7,8,8,9,9,9,10,10,10],
	[5,3,6,9,2,6,5,8,1,10,6,9,3,4,7,5,9,1,3,7,2,6,10],ones(Int,23)) 
    #rectangular matrix
    R = sparse([1,1,1,2,2,2,3,3,4,4,5,5,6,6,6,7,8,8,8,9,9,10],
    [1,2,4,2,3,5,1,3,2,4,3,5,2,3,6,4,1,2,6,3,4,3],ones(Int,22))

    p1,p2,p3 = etree(D),etree(C),etree(R)
    porder1,porder2,porder3 = postorder(p1),postorder(p2),postorder(p3)

@testset "etree" begin  
    @test p1 == [4,3,4,7,6,8,8,9,0]
    @test p2 == [3,6,4,6,8,7,9,9,10,0]
    @test p3 == [2,3,4,5,6,0]
end
@testset "postorder" begin
    @test porder1 == [5,6,1,2,3,4,7,8,9]
    @test porder2 == [2,1,3,4,6,7,5,8,9,10]
    @test porder3 == [1,2,3,4,5,6]
end
    nrc1 = rowcolcount!(D,p1,porder1)
    rc1 = rowcolcount!(D,p1,porder1,true)

    nrc2,cc2 = rowcolcount!(C,p2,porder2)
    rc2,cc2 = rowcolcount!(C,p2,porder2,true)

    nrc3,cc3 = rowcolcount!(R,p3,porder3)
    rc3,cc3 = rowcolcount!(R,p3,porder3,true)
@testset "rowcolcount!" begin
    @test nrc1 == [3,4,4,3,3,2,3,2,1]
    @test rc1 == [3,4,4,3,3,0,3,2,0]
    @test p1 == [4,3,4,7,8,-5,8,0,-8]

    @test nrc2 == [4,3,6,5,3,4,3,2,2,1]
    @test rc2 == [4,3,6,0,3,4,0,0,2,0]
    @test cc2 == [2,2,3,2,3,3,2,2,2,1]
    @test p2 == [3,6,6,-3,9,9,-6,-5,0,-9]

    @test nrc3 == [5,5,4,3,2,1]
    @test rc3 == [5,5,0,0,0,0]
    @test cc3 == [3,5,7,7,6,5]
    @test p3 == [2,0,-2,-2,-2,-2]
end
    np1 = supernodecount(p1)
    np2 = supernodecount(p2)
    np3 = supernodecount(p3)
@testset "supernodecount" begin
    @test np1 == [1,1,1,1,2,0,1,2,0]
    @test np2 == [1,1,2,0,2,2,0,0,2,0]
    @test np3 == [1,5,0,0,0,0]
end
    idxr1 = idx_r(D,p1,porder1)
@testset "idx_r" begin
    @test idxr1 == [Set([1,4,9]),Set([2,3,4,7]),Set([3,4,7,9]),Set([4,7,9]),Set([5,6,8]),
                    Set([]),Set([7,8,9]),Set([8,9]),Set([])]

    @test_throws AssertionError idx_r(C,p2,porder2)
    @test_throws AssertionError idx_r(R,p3,porder3)
end
@testset "symbolicfact" begin
    s1 = symbolicfact(D)
    s2 = symbolicfact(C)
    s3 = symbolicfact(R)

    for i = 1:length(s1)
        @test s1[i].father == p1[i]
        @test s1[i].rc == rc1[i]
        @test s1[i].cc == -1
        @test s1[i].idxr == idxr1[i]
        @test s1[i].idxc == Set([])
        @test s1[i].np == np1[i]
    end

    for i = 1:length(s2)
        @test s2[i].father == p2[i]
        @test s2[i].rc == rc2[i]
        @test s2[i].cc == cc2[i]
        @test s2[i].idxr == Set([])
        @test s2[i].idxc == Set([])
        @test s2[i].np == np2[i]
    end

    for i = 1:length(s3)
        @test s3[i].father == p3[i]
        @test s3[i].rc == rc3[i]
        @test s3[i].cc == cc3[i]
        @test s3[i].idxr == Set([])
        @test s3[i].idxc == Set([])
        @test s3[i].np == np3[i]
    end
end