
@enum MatOperation begin
  MATOP_SET_VALUES=0
  MATOP_GET_ROW=1
  MATOP_RESTORE_ROW=2
  MATOP_MULT=3
  MATOP_MULT_ADD=4
  MATOP_MULT_TRANSPOSE=5
  MATOP_MULT_TRANSPOSE_ADD=6
  MATOP_SOLVE=7
  MATOP_SOLVE_ADD=8
  MATOP_SOLVE_TRANSPOSE=9
  MATOP_SOLVE_TRANSPOSE_ADD=10
  MATOP_LUFACTOR=11
  MATOP_CHOLESKYFACTOR=12
  MATOP_SOR=13
  MATOP_TRANSPOSE=14
  MATOP_GETINFO=15
  MATOP_EQUAL=16
  MATOP_GET_DIAGONAL=17
  MATOP_DIAGONAL_SCALE=18
  MATOP_NORM=19
  MATOP_ASSEMBLY_BEGIN=20
  MATOP_ASSEMBLY_END=21
  MATOP_SET_OPTION=22
  MATOP_ZERO_ENTRIES=23
  MATOP_ZERO_ROWS=24
  MATOP_LUFACTOR_SYMBOLIC=25
  MATOP_LUFACTOR_NUMERIC=26
  MATOP_CHOLESKY_FACTOR_SYMBOLIC=27
  MATOP_CHOLESKY_FACTOR_NUMERIC=28
  MATOP_SETUP_PREALLOCATION=29
  MATOP_ILUFACTOR_SYMBOLIC=30
  MATOP_ICCFACTOR_SYMBOLIC=31
  MATOP_GET_DIAGONAL_BLOCK=32
  MATOP_FREE_INTER_STRUCT=33
  MATOP_DUPLICATE=34
  MATOP_FORWARD_SOLVE=35
  MATOP_BACKWARD_SOLVE=36
  MATOP_ILUFACTOR=37
  MATOP_ICCFACTOR=38
  MATOP_AXPY=39
  MATOP_CREATE_SUBMATRICES=40
  MATOP_INCREASE_OVERLAP=41
  MATOP_GET_VALUES=42
  MATOP_COPY=43
  MATOP_GET_ROW_MAX=44
  MATOP_SCALE=45
  MATOP_SHIFT=46
  MATOP_DIAGONAL_SET=47
  MATOP_ZERO_ROWS_COLUMNS=48
  MATOP_SET_RANDOM=49
  MATOP_GET_ROW_IJ=50
  MATOP_RESTORE_ROW_IJ=51
  MATOP_GET_COLUMN_IJ=52
  MATOP_RESTORE_COLUMN_IJ=53
  MATOP_FDCOLORING_CREATE=54
  MATOP_COLORING_PATCH=55
  MATOP_SET_UNFACTORED=56
  MATOP_PERMUTE=57
  MATOP_SET_VALUES_BLOCKED=58
  MATOP_CREATE_SUBMATRIX=59
  MATOP_DESTROY=60
  MATOP_VIEW=61
  MATOP_CONVERT_FROM=62
  MATOP_MATMAT_MULT=63
  MATOP_MATMAT_MULT_SYMBOLIC=64
  MATOP_MATMAT_MULT_NUMERIC=65
  MATOP_SET_LOCAL_TO_GLOBAL_MAP=66
  MATOP_SET_VALUES_LOCAL=67
  MATOP_ZERO_ROWS_LOCAL=68
  MATOP_GET_ROW_MAX_ABS=69
  MATOP_GET_ROW_MIN_ABS=70
  MATOP_CONVERT=71
  MATOP_SET_COLORING=72
  # /* MATOP_PLACEHOLDER_73=73 */
  MATOP_SET_VALUES_ADIFOR=74
  MATOP_FD_COLORING_APPLY=75
  MATOP_SET_FROM_OPTIONS=76
  MATOP_MULT_CONSTRAINED=77
  MATOP_MULT_TRANSPOSE_CONSTRAIN=78
  MATOP_FIND_ZERO_DIAGONALS=79
  MATOP_MULT_MULTIPLE=80
  MATOP_SOLVE_MULTIPLE=81
  MATOP_GET_INERTIA=82
  MATOP_LOAD=83
  MATOP_IS_SYMMETRIC=84
  MATOP_IS_HERMITIAN=85
  MATOP_IS_STRUCTURALLY_SYMMETRIC=86
  MATOP_SET_VALUES_BLOCKEDLOCAL=87
  MATOP_CREATE_VECS=88
  MATOP_MAT_MULT=89
  MATOP_MAT_MULT_SYMBOLIC=90
  MATOP_MAT_MULT_NUMERIC=91
  MATOP_PTAP=92
  MATOP_PTAP_SYMBOLIC=93
  MATOP_PTAP_NUMERIC=94
  MATOP_MAT_TRANSPOSE_MULT=95
  MATOP_MAT_TRANSPOSE_MULT_SYMBO=96
  MATOP_MAT_TRANSPOSE_MULT_NUMER=97
  # /* MATOP_PLACEHOLDER_98=98 */
  MATOP_PRODUCTSETFROMOPTIONS=99
  MATOP_PRODUCTSYMBOLIC=100
  MATOP_PRODUCTNUMERIC=101
  MATOP_CONJUGATE=102
  # /* MATOP_PLACEHOLDER_103=103 */
  MATOP_SET_VALUES_ROW=104
  MATOP_REAL_PART=105
  MATOP_IMAGINARY_PART=106
  MATOP_GET_ROW_UPPER_TRIANGULAR=107
  MATOP_RESTORE_ROW_UPPER_TRIANG=108
  MATOP_MAT_SOLVE=109
  MATOP_MAT_SOLVE_TRANSPOSE=110
  MATOP_GET_ROW_MIN=111
  MATOP_GET_COLUMN_VECTOR=112
  MATOP_MISSING_DIAGONAL=113
  MATOP_GET_SEQ_NONZERO_STRUCTUR=114
  MATOP_CREATE=115
  MATOP_GET_GHOSTS=116
  MATOP_GET_LOCAL_SUB_MATRIX=117
  MATOP_RESTORE_LOCALSUB_MATRIX=118
  MATOP_MULT_DIAGONAL_BLOCK=119
  MATOP_HERMITIAN_TRANSPOSE=120
  MATOP_MULT_HERMITIAN_TRANSPOSE=121
  MATOP_MULT_HERMITIAN_TRANS_ADD=122
  MATOP_GET_MULTI_PROC_BLOCK=123
  MATOP_FIND_NONZERO_ROWS=124
  MATOP_GET_COLUMN_NORMS=125
  MATOP_INVERT_BLOCK_DIAGONAL=126
  # /* MATOP_PLACEHOLDER_127=127 */
  MATOP_CREATE_SUB_MATRICES_MPI=128
  MATOP_SET_VALUES_BATCH=129
  MATOP_TRANSPOSE_MAT_MULT=130
  MATOP_TRANSPOSE_MAT_MULT_SYMBO=131
  MATOP_TRANSPOSE_MAT_MULT_NUMER=132
  MATOP_TRANSPOSE_COLORING_CREAT=133
  MATOP_TRANS_COLORING_APPLY_SPT=134
  MATOP_TRANS_COLORING_APPLY_DEN=135
  MATOP_RART=136
  MATOP_RART_SYMBOLIC=137
  MATOP_RART_NUMERIC=138
  MATOP_SET_BLOCK_SIZES=139
  MATOP_AYPX=140
  MATOP_RESIDUAL=141
  MATOP_FDCOLORING_SETUP=142
  MATOP_MPICONCATENATESEQ=144
  MATOP_DESTROYSUBMATRICES=145
  MATOP_TRANSPOSE_SOLVE=146
  MATOP_GET_VALUES_LOCAL=147
end

macro wrapper_bis(fn,rt,argts,args,url)
  hn = Symbol("$(fn.value)_handle")
  sargs = "$(args)"
  if length(args.args) == 1
    sargs = sargs[1:end-2]*")"
  end
  if isempty(rstrip(url))
    str = """
        $(fn.value)$(sargs)
    """
  else
    str = """
        $(fn.value)$(sargs)

    See [PETSc manual]($url).
    """
  end
  expr = quote
    const $hn = Ref(C_NULL)
    push!(GridapPETSc._PRELOADS,($hn,$fn))
    @doc $str
    function $(fn.value)($(args.args...))
      _handle = Libdl.dlsym(GridapPETSc.PETSC.libpetsc_handle[],$fn;throw_error=false)
      $(hn)[]=_handle
      @check $(hn)[] != C_NULL "Missing symbol. Re-configure and compile PETSc."
      ccall($(hn)[],$rt,$argts,$(args.args...))
    end
  end
  esc(expr)
end

# MatShell
@wrapper_bis(:MatCreateShell,PetscErrorCode,(MPI.Comm,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Cvoid},Ptr{Mat}),(comm,m,n,M,N,ctx,mat),"https://petsc.org/main/docs/manualpages/Mat/MatCreateShell/")
@wrapper_bis(:MatShellSetOperation,PetscErrorCode,(Mat,MatOperation,Ptr{Cvoid}),(mat,matop,g),"https://petsc.org/release/docs/manualpages/Mat/MatShellSetOperation.html")
@wrapper_bis(:MatShellGetContext,PetscErrorCode,(Mat,Ptr{Ptr{Cvoid}}),(mat,ctx),"https://petsc.org/release/docs/manualpages/Mat/MatShellGetContext.html")
# PCMG
@wrapper_bis(:PCMGSetLevels,PetscErrorCode,(PC,PetscInt,Ptr{MPI.Comm}),(pc,levels,comms),"https://petsc.org/main/docs/manualpages/PC/PCMGSetLevels/")
@wrapper_bis(:PCMGSetResidual,PetscErrorCode,(PC,PetscInt,Ptr{Cvoid},Mat),(pc,l,residual,mat),"https://petsc.org/main/docs/manualpages/PC/PCMGSetResidual/")
@wrapper_bis(:PCMGSetInterpolation,PetscErrorCode,(PC,PetscInt,Mat),(pc,l,mat),"https://petsc.org/release/docs/manualpages/PC/PCMGSetInterpolation.html")
@wrapper_bis(:PCMGSetRestriction,PetscErrorCode,(PC,PetscInt,Mat),(pc,l,mat),"https://petsc.org/release/docs/manualpages/PC/PCMGSetRestriction.html")
@wrapper_bis(:PCMGGetSmoother, PetscErrorCode, (PC,PetscInt,Ptr{KSP}), (pc,l,ksp), "https://petsc.org/release/docs/manualpages/PC/PCMGGetSmoother.html")
@wrapper_bis(:PCMGGetCoarseSolve, PetscErrorCode, (PC,Ptr{KSP}), (pc,ksp), "https://petsc.org/release/docs/manualpages/PC/PCMGGetCoarseSolve.html")


# DMPlex

"""
Julia alias for the `DM` C type.

See [PETSc manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DM/DM.html).
"""
struct DM
  ptr::Ptr{Cvoid}
end
DM() = Vec(Ptr{Cvoid}())
Base.convert(::Type{DM},p::Ptr{Cvoid}) = DM(p)
Base.unsafe_convert(::Type{Ptr{Cvoid}},v::DM) = v.ptr

# PetscErrorCode DMPlexCreate(MPI_Comm comm, DM *mesh)
@wrapper_bis(:DMPlexCreate, PetscErrorCode, (MPI.Comm,Ptr{DM}), (comm,mesh), "https://petsc.org/release/docs/manualpages/DMPLEX/DMPlexCreate.html")

# PetscErrorCode  DMDestroy(DM *dm)
@wrapper_bis(:DMDestroy, PetscErrorCode, (Ptr{DM},), (dm,), "https://petsc.org/release/docs/manualpages/DM/DMDestroy.html")

@wrapper_bis(:DMView,PetscErrorCode,(DM,GridapPETSc.PETSC.PetscViewer),(dm,viewer),"https://petsc.org/release/docs/manualpages/DM/DMView.html")

# PetscErrorCode DMPlexSetChart(DM dm, PetscInt pStart, PetscInt pEnd)
@wrapper_bis(:DMPlexSetChart, PetscErrorCode, (DM,PetscInt,PetscInt), (dm,pStart,pEnd), "https://petsc.org/release/docs/manualpages/DMPLEX/DMPlexSetChart.html")

@wrapper_bis(:DMPlexGetChart,PetscErrorCode,(DM,Ptr{PetscInt},Ptr{PetscInt}),(dm,pStart,pEnd),"https://petsc.org/release/docs/manualpages/DMPLEX/DMPlexGetChart.html")

# PetscErrorCode DMPlexGetHeightStratum(DM dm, PetscInt stratumValue, PetscInt *start, PetscInt *end)
@wrapper_bis(:DMPlexGetHeightStratum,PetscErrorCode,(DM,PetscInt,Ptr{PetscInt},Ptr{PetscInt}),(dm,stratumValue,s,e),"https://petsc.org/release/docs/manualpages/DMPLEX/DMPlexGetHeightStratum.html")

# PetscErrorCode DMPlexSetConeSize(DM dm, PetscInt p, PetscInt size)
@wrapper_bis(:DMPlexSetConeSize, PetscErrorCode, (DM,PetscInt,PetscInt), (dm,p,size),  "https://petsc.org/release/docs/manualpages/DMPLEX/DMPlexSetConeSize.html")

# PetscErrorCode DMSetUp(DM dm)
@wrapper_bis(:DMSetUp, PetscErrorCode, (DM,), (dm,), "https://petsc.org/release/docs/manualpages/DM/DMSetUp.html")

# PetscErrorCode DMPlexSetCone(DM dm, PetscInt p, const PetscInt cone[])
@wrapper_bis(:DMPlexSetCone, PetscErrorCode, (DM,PetscInt,Ptr{PetscInt}), (dm,p,cone), "https://petsc.org/release/docs/manualpages/DMPLEX/DMPlexSetCone.html")

# PetscErrorCode DMPlexSymmetrize(DM dm)
@wrapper_bis(:DMPlexSymmetrize, PetscErrorCode, (DM,), (dm,), "https://petsc.org/release/docs/manualpages/DMPLEX/DMPlexSymmetrize.html")

# PetscErrorCode DMPlexStratify(DM dm)
@wrapper_bis(:DMPlexStratify, PetscErrorCode, (DM,), (dm,), "https://petsc.org/release/docs/manualpages/DMPLEX/DMPlexStratify.html")

@wrapper_bis(:DMPlexCheck, PetscErrorCode, (DM,), (dm,), "https://petsc.org/release/docs/manualpages/DMPLEX/DMPlexCheck.html")


"""
Julia alias for the `PetscSection` C type.

See [PETSc manual](https://petsc.org/release/docs/manualpages/PetscSection/PetscSection.html).
"""
struct PetscSection
  ptr::Ptr{Cvoid}
end
PetscSection() = Vec(Ptr{Cvoid}())
Base.convert(::Type{PetscSection},p::Ptr{Cvoid}) = PetscSection(p)
Base.unsafe_convert(::Type{Ptr{Cvoid}},v::PetscSection) = v.ptr

# PetscSectionCreate(MPI_Comm,PetscSection *);
@wrapper_bis(:PetscSectionCreate, PetscErrorCode, (MPI.Comm,Ptr{PetscSection}), (comm,sec), "https://petsc.org/release/docs/manualpages/PetscSection/PetscSectionCreate.html")

# PetscSectionSetChart(PetscSection,low,high);
@wrapper_bis(:PetscSectionSetChart, PetscErrorCode, (PetscSection,PetscInt,PetscInt), (sec,low,high), "https://petsc.org/release/docs/manualpages/PetscSection/PetscSectionSetChart.html")

# PetscSectionSetDof(PetscSection,point,numdof);
@wrapper_bis(:PetscSectionSetDof, PetscErrorCode, (PetscSection,PetscInt,PetscInt), (sec,point,numdof), "https://petsc.org/release/docs/manualpages/PetscSection/PetscSectionSetDof.html")

# PetscSectionSetUp(PetscSection);
@wrapper_bis(:PetscSectionSetUp, PetscErrorCode, (PetscSection,), (sec,), "https://petsc.org/release/docs/manualpages/PetscSection/PetscSectionSetUp.html")

# PetscSectionSetNumFields(PetscSection, numFields);
@wrapper_bis(:PetscSectionSetNumFields, PetscErrorCode, (PetscSection,PetscInt), (sec,numFields), "https://petsc.org/release/docs/manualpages/PetscSection/PetscSectionSetNumFields.html")

# PetscSectionGetOffset(PetscSection,point,PetscInt *);
@wrapper_bis(:PetscSectionGetOffset, PetscErrorCode, (PetscSection,PetscInt,Ptr{PetscInt}), (sec,point,off), "https://petsc.org/release/docs/manualpages/PetscSection/PetscSectionGetOffset.html")

# PetscErrorCode PetscSectionSetConstraintDof(PetscSection s, PetscInt point, PetscInt numDof)
@wrapper_bis(:PetscSectionSetConstraintDof, PetscErrorCode, (PetscSection,PetscInt,PetscInt), (s,point,numDof), "https://petsc.org/release/docs/manualpages/PetscSection/PetscSectionSetConstraintDof.html")

# PetscSectionDestroy(PetscSection);
@wrapper_bis(:PetscSectionDestroy, PetscErrorCode, (Ptr{PetscSection},), (sec,), "https://petsc.org/release/docs/manualpages/PetscSection/PetscSectionDestroy.html")
