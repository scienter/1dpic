EXEC = show
CC = mpicc 
OBJS = main.o boundary.o loadPlasma_random.o parameterSetting.o findparam.o saveFile.o loadLaser.o fieldSolve.o interpolation.o particlePush.o updateCurrent.o removeEdge.o movingDomain.o fieldShareX.o filter.o boostShot.o probe.o rearrangeParticles.o clean.o particleShareX.o dumpData.o loadPlasma_crystal.o fieldIonization.o
INCL = constants.h laser.h mesh.h particle.h plasma.h
LIBS = -lm
$(EXEC):$(OBJS)
	$(CC) $(OBJS) $(LIBS) -o $(EXEC)
$(OBJS):$(INCL)
clean:
	@rm -f *.o *~ $(EXEC)
