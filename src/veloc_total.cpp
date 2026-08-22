#include <petscksp.h>
#include <petscmath.h>

PetscErrorCode build_veloc(int dimensions);

PetscErrorCode solve_veloc(int dimensions);

extern Vec Veloc;
extern Vec Veloc_fut;
extern Vec Veloc_step1;
extern Vec Veloc_step2;

extern PetscInt WITH_NON_LINEAR;

extern int tcont;
extern double visc_MAX;
extern double visc_MIN;

extern PetscReal Xi_min;

extern Vec Veloc_0, Veloc_0_copy;
int step_nl = 0;
extern PetscBool init_winkler;
extern PetscReal c_winkler;
extern PetscErrorCode calc_kinematic_winkler();

PetscErrorCode veloc_total(int dimensions)
{
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	// winkler to modify basal b.c. (Veloc_0 and Veloc_fut) with the restauration velocity

	if (dimensions==2 && c_winkler > 0 && init_winkler == PETSC_TRUE) {
			VecCopy(Veloc_0,Veloc_0_copy);
			ierr = calc_kinematic_winkler(); CHKERRQ(ierr);
	}

	PetscPrintf(PETSC_COMM_WORLD,"\nSolving Stokes system:\n");

	PetscInt n;

	VecCopy(Veloc_fut,Veloc);

	if (WITH_NON_LINEAR==0){
		ierr = build_veloc(dimensions);CHKERRQ(ierr);

		ierr = solve_veloc(dimensions);CHKERRQ(ierr);
	}
	else {

		VecCopy(Veloc_fut,Veloc_step1);

		VecGetSize(Veloc_fut,&n);


		PetscReal VM1,VM2,sig1,sig2,vivi;

		PetscReal Xi=50000.0;

		for (int step_nl=0; step_nl<700 && Xi>Xi_min; step_nl++){

			
			ierr = build_veloc(dimensions);CHKERRQ(ierr);

			ierr = solve_veloc(dimensions);CHKERRQ(ierr);

			VecCopy(Veloc_fut,Veloc_step2);


			VecSum(Veloc_step1,&VM1);
			VM1/=n;
			VecShift(Veloc_step1,-VM1);


			VecSum(Veloc_step2,&VM2);
			VM2/=n;
			VecShift(Veloc_step2,-VM2);



			VecDot(Veloc_step1,Veloc_step1,&sig1);
			VecDot(Veloc_step2,Veloc_step2,&sig2);

			VecDot(Veloc_step1,Veloc_step2,&vivi);
			//vivi/=n;

			Xi = 1.0 - vivi/PetscSqrtReal(sig1*sig2);

			PetscPrintf(PETSC_COMM_WORLD,"      Xi = %lg %d\n\n",Xi,step_nl);

			VecCopy(Veloc_fut,Veloc_step1);

		}
	}

	VecCopy(Veloc_0_copy,Veloc_0);

	PetscFunctionReturn(0);

}
