#include<iostream>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/LU>
#include"petscksp.h"
#include<vector>

using namespace Eigen;

MatrixXd solve_linear(
    const SparseMatrix<double> A,
    const MatrixXd rhs,
	const int linear_solver_type,
	const double tolerance)
{
	int A_rows = A.rows();
	int rhs_cols = rhs.cols();
    MatrixXd X(A_rows, rhs_cols);
    X.setZero();

	assert(A_rows == A.cols());
	assert(A_rows == rhs.rows());
    if (linear_solver_type == 0) {
		// Sparse LU
        SparseLU <SparseMatrix<double>, COLAMDOrdering< int > > slusolver;
        slusolver.analyzePattern(A);
        slusolver.factorize(A);
        if (slusolver.info() != 0)
            std::cout<<"Factorization failed. Error: "<<slusolver.info()<<std::endl;
        X = slusolver.solve(rhs);
    } else if (linear_solver_type == 1) {
		// Dense LU full pivoting
		MatrixXd matAdense(A_rows, A_rows);
		MatrixXd eye(A_rows, A_rows);
		eye.setIdentity();
		matAdense = A * eye;

		//double offset = 0;//0.00001;
		//matAdense = matAdense + eye * offset;
        // Full Pivoting LU Factorization
        X = matAdense.fullPivLu().solve(rhs);
    } else if (linear_solver_type == 2) {
		// BiCGSTAB
        BiCGSTAB<SparseMatrix<double> > itsolver;

		itsolver.setTolerance(tolerance);
        itsolver.compute(A);
        if (itsolver.info() == 0) std::cout<<"Iterative Factorization success"<<std::endl;
        else std::cout<<"Factorization failed. Error: "<<itsolver.info()<<std::endl;

        X = itsolver.solve(rhs);
        std::cout << "#iterations:     " << itsolver.iterations() << std::endl;
        std::cout << "estimated error: " << itsolver.error()      << std::endl;

    } else if (linear_solver_type == 3) {

		Vec         x_petsc, b_petsc;
		Mat         A_petsc;
		KSP         ksp;
		PC          pc;
		PetscInt    m = A.rows(), n = A.cols();
		PetscInt    r, c;
		PetscScalar v;
		PetscInt    *indices;
		PetscScalar *values;

	//  MatCreate(PETSC_COMM_SELF,&A);
	//  MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m,n);
	//  MatSetFromOptions(A);
	//  MatSetUp(A);
		MatCreateSeqAIJ(PETSC_COMM_SELF, m, n, 9, 0, &A_petsc);

		for (int k = 0; k < A.outerSize(); ++k)
		{
			// Iterate over inside
			for (SparseMatrix<double>::InnerIterator it(A,k); it; ++it)
			{
				r = it.row();
				c = it.col();
				v = it.value();
				MatSetValues(A_petsc, 1, &r, 1, &c, &v, INSERT_VALUES);
			}
		}
		MatAssemblyBegin(A_petsc,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(A_petsc,MAT_FINAL_ASSEMBLY);


		KSPCreate(PETSC_COMM_SELF,&ksp);

		KSPSetOperators(ksp,A_petsc,A_petsc);
		KSPGetPC(ksp,&pc);
		PCSetType(pc,PCNONE);
		//PCSetType(pc,PCSOR);
		//PCSetType(pc,PCJACOBI);
		//PCSetType(pc,PCILU);

		KSPSetType(ksp, KSPGMRES);
		KSPGMRESSetRestart(ksp, 300);
	//  KSPSetTolerances(ksp,1.e-1,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
		KSPSetTolerances(ksp,tolerance,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
	//  KSPSetMonitor(ksp,KSPDefaultMonitor,PETSC_NULL,PETSC_NULL);
		PetscViewerAndFormat *vf;
		PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,&vf);
		KSPMonitorSet(ksp,(PetscErrorCode (*)(KSP,PetscInt,PetscReal,void*))KSPMonitorTrueResidualNorm,vf,(PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy);
		KSPSetFromOptions(ksp);


		PetscMalloc1(rhs.rows() * sizeof(PetscScalar), &values);
		PetscMalloc1(rhs.rows() * sizeof(PetscInt), &indices);
		VecCreate(PETSC_COMM_SELF, &b_petsc);
		VecSetSizes(b_petsc,PETSC_DECIDE, rhs.rows());
		VecSetFromOptions(b_petsc);
		VecDuplicate(b_petsc, &x_petsc);
		for (int bcol = 0; bcol < rhs.cols(); bcol++)
		{

			for (int brow = 0; brow < rhs.rows(); brow++)
			{
				v = rhs(brow, bcol);
				VecSetValue(b_petsc, brow, v, INSERT_VALUES);
			}

			KSPSolve(ksp,b_petsc,x_petsc);


			for (PetscInt brow = 0; brow < rhs.rows(); brow++)
			{
				values[brow] = -1.0;
				indices[brow] = brow;
			}
			VecGetValues(x_petsc, rhs.rows(), indices, values);

			for (int brow = 0; brow < rhs.rows(); brow++)
			{
				X(brow, bcol) = values[brow];
			}

		}
		PetscFree(values);
		PetscFree(indices);
		VecDestroy(&x_petsc);
		VecDestroy(&b_petsc);
		MatDestroy(&A_petsc);

		KSPDestroy(&ksp);
	}

    return X;
}


VectorXd solve_linear(
    const SparseMatrix<double> A,
    const VectorXd rhs,
	const int linear_solver_type,
	const double tolerance)
{
	int A_rows = A.rows();
	int A_cols = A.cols();
	int rhs_rows = rhs.rows();
	assert(A_rows==A_cols);
	assert(A_rows==rhs_rows);
    VectorXd X(rhs_rows);
    X.setZero();

    if (linear_solver_type == 0) {
		// Sparse LU
        SparseLU <SparseMatrix<double>, COLAMDOrdering< int > > slusolver;
        slusolver.analyzePattern(A);
        slusolver.factorize(A);
        if (slusolver.info() != 0)
            std::cout<<"Factorization failed. Error: "<<slusolver.info()<<std::endl;
        X = slusolver.solve(rhs);
    } else if (linear_solver_type == 1) {
		// Dense LU full pivoting
		MatrixXd matAdense(A_rows, A_rows);
		MatrixXd eye(A_rows, A_rows);
		eye.setIdentity();
		matAdense = A * eye;

		//double offset = 0;//0.00001;
		//matAdense = matAdense + eye * offset;
        // Full Pivoting LU Factorization
        X = matAdense.fullPivLu().solve(rhs);
    } else if (linear_solver_type == 2) {
		// BiCGSTAB
        BiCGSTAB<SparseMatrix<double> > itsolver;
        itsolver.compute(A);
        if (itsolver.info() == 0)
            std::cout<<"Iterative Factorization success"<<std::endl;
        else
            std::cout<<"Factorization failed. Error: "<<itsolver.info()<<std::endl;
        std::cout << "#iterations:     " << itsolver.iterations() << std::endl;
        std::cout << "estimated error: " << itsolver.error()      << std::endl;
        X = itsolver.solve(rhs);
    } else if (linear_solver_type == 3) {
		Vec         x_petsc, rhs_petsc;
		Mat         A_petsc;
		KSP         ksp;
		PC          pc;
		PetscInt    m = A.rows(), n = A.cols();
		PetscInt    r, c;
		PetscScalar v;
		PetscInt    *indices;
		PetscScalar *values;

	//  MatCreate(PETSC_COMM_SELF,&A);
	//  MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m,n);
	//  MatSetFromOptions(A);
	//  MatSetUp(A);
		MatCreateSeqAIJ(PETSC_COMM_SELF, m, n, 9, 0, &A_petsc);

		for (int k = 0; k < A.outerSize(); ++k)
		{
			// Iterate over inside
			for (SparseMatrix<double>::InnerIterator it(A,k); it; ++it)
			{
				r = it.row();
				c = it.col();
				v = it.value();
				MatSetValues(A_petsc, 1, &r, 1, &c, &v, INSERT_VALUES);
			}
		}
		MatAssemblyBegin(A_petsc,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(A_petsc,MAT_FINAL_ASSEMBLY);


		KSPCreate(PETSC_COMM_SELF,&ksp);

		KSPSetOperators(ksp,A_petsc,A_petsc);
		KSPGetPC(ksp,&pc);
		PCSetType(pc,PCNONE);
	//  PCSetType(pc,PCJACOBI);
	//  PCSetType(pc,PCILU);

		KSPSetType(ksp, KSPGMRES);
		KSPGMRESSetRestart(ksp, 300);
	//  KSPSetTolerances(ksp,1.e-1,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
		KSPSetTolerances(ksp,tolerance,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
	//  KSPSetMonitor(ksp,KSPDefaultMonitor,PETSC_NULL,PETSC_NULL);
		PetscViewerAndFormat *vf;
		PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,&vf);
		KSPMonitorSet(ksp,(PetscErrorCode (*)(KSP,PetscInt,PetscReal,void*))KSPMonitorTrueResidualNorm,vf,(PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy);
		KSPSetFromOptions(ksp);


		PetscMalloc1(rhs.rows() * sizeof(PetscScalar), &values);
		PetscMalloc1(rhs.rows() * sizeof(PetscInt), &indices);
		VecCreate(PETSC_COMM_SELF, &rhs_petsc);
		VecSetSizes(rhs_petsc,PETSC_DECIDE, rhs.rows());
		VecSetFromOptions(rhs_petsc);
		VecDuplicate(rhs_petsc, &x_petsc);
		for (int bcol = 0; bcol < rhs.cols(); bcol++)
		{

			for (int brow = 0; brow < rhs.rows(); brow++)
			{
				v = rhs(brow, bcol);
				VecSetValue(rhs_petsc, brow, v, INSERT_VALUES);
			}

			KSPSolve(ksp,rhs_petsc,x_petsc);


			for (PetscInt brow = 0; brow < rhs.rows(); brow++)
			{
				values[brow] = -1.0;
				indices[brow] = brow;
			}
			VecGetValues(x_petsc, rhs.rows(), indices, values);

			for (int brow = 0; brow < rhs.rows(); brow++)
			{
				X(brow, bcol) = values[brow];
			}

		}
		PetscFree(values);
		PetscFree(indices);
		VecDestroy(&x_petsc);
		VecDestroy(&rhs_petsc);
		MatDestroy(&A_petsc);

		KSPDestroy(&ksp);

	}

    return X;
}

//VectorXd solve_linear(
//    const MatrixXd A,
//    const VectorXd rhs,
//	const int linear_solver_type,
//	const double tolerance)
//{
//	int A_rows = A.rows();
//	int A_cols = A.cols();
//	int rhs_rows = rhs.rows();
//	assert(A_rows==A_cols);
//	assert(A_rows==rhs_rows);
//    VectorXd X(rhs_rows);
//    X.setZero();
//
//    if (linear_solver_type == 0) {
//        X = A.fullPivLu().solve(rhs);
//    } else {
//        abort();
//	}
//
//    return X;
//}

MatrixXd solve_dense_linear(
    const MatrixXd A,
    const MatrixXd rhs,
	const int linear_solver_type,
	const double tolerance)
{
	int A_rows = A.rows();
	int A_cols = A.cols();
	int rhs_rows = rhs.rows();
	int rhs_cols = rhs.cols();
	assert(A_rows==A_cols);
	assert(A_rows==rhs_rows);
    MatrixXd X(rhs_rows,rhs_cols);
    X.setZero();

    if (linear_solver_type == 0) {
        FullPivLU<MatrixXd> lu(A);
		lu.setThreshold(1.0e-13);
        const int mrank = lu.rank();
		std::cout<<mrank<<"    "<<A_rows<<std::endl;
	std::cout<<A<<std::endl<<std::endl<<std::endl;
	std::cout<<rhs<<std::endl<<std::endl<<std::endl;
		if (mrank < A_rows) {
            std::cout<<"Dense LU failed. "<<std::endl;
			abort();
		} else {
			X = lu.solve(rhs);
	std::cout<<X<<std::endl<<std::endl<<std::endl;
		}
    } else {
        abort();
	}

    return X;
}
