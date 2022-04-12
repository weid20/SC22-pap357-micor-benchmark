int get_logp(unsigned int value)
{
	int ret = 0;
	int tmp = 0x1;
	while(tmp < value)
	{
		ret++;
		tmp <<= 1;
	}
	return ret;
}

int
get_bwr(int value, int bit)
{
	int i;
	int ret = 0;
	int tmp = 0x1;
	
	for(i = 0; i < bit; i++)
	{
		if(value & (0x1<<(bit - i - 1)))
			ret += tmp;
		tmp <<= 1;
	}
	return ret;
}

#undef FUNCNAME
#define FUNCNAME MPIC_Gather_positive
#undef FCNAME
#define FCNAME MPIU_QUOTE(FUNCNAME)
int MPIC_Gather_positive(
		const void *sendbuf,
		int sendbytes,
		void *recvbuf,
		int recvbytes,
		int root,
		MPID_Comm *comm_ptr,
		int need_shuffle,
		int *errflag)
{
    MPI_Status status;
    int        rank, comm_size;
    int mask, src, dst;
    int        mpi_errno = MPI_SUCCESS;
	char *send_start, *recv_start;
    MPI_Comm comm;
	int logp;
	MPI_Status sta;
    MPIU_CHKLMEM_DECL(1);
	

	int nbytes;
	int j = 0;
    
    comm = comm_ptr->handle;
    comm_size = comm_ptr->local_size;
    rank = comm_ptr->rank;
	logp = get_logp(comm_size);
	mask = 0x1;

	dst = get_bwr(rank, logp);
	
	if(need_shuffle == 1 && rank != dst)
	{
		send_start = (char *)recvbuf + rank * recvbytes;
		MPIC_Sendrecv(sendbuf, sendbytes, MPI_BYTE, dst, MPIR_GATHER_TAG, send_start, recvbytes, MPI_BYTE, dst, MPIR_GATHER_TAG, comm_ptr, &sta, errflag);
	}
	else
	{
		send_start = sendbuf;
	}

	nbytes = sendbytes;

	/* Use binomial-tree algorithm */
	while(mask < comm_size)
	{
		if((rank & (mask - 1)) == 0)
		{
			dst = rank ^ mask;
			recv_start = send_start + nbytes;
			if(rank < dst)
			{
				MPIC_Recv(recv_start, nbytes, MPI_BYTE, dst, MPIR_GATHER_TAG, comm_ptr, &sta, errflag);
			}
			else
			{
				MPIC_Send(send_start, nbytes, MPI_BYTE, dst, MPIR_GATHER_TAG, comm_ptr, errflag);
			}
		}
		mask <<= 1;
		nbytes <<= 1;
		j++;
	}

fn_exit:
    return mpi_errno;
 fn_fail:
    goto fn_exit;

}


#undef FUNCNAME
#define FUNCNAME MPIC_Gather_negative
#undef FCNAME
#define FCNAME MPIU_QUOTE(FUNCNAME)
int MPIC_Gather_negative(
		const void *sendbuf,
		int sendbytes,
		void *recvbuf,
		int recvbytes,
		int root,
		MPID_Comm *comm_ptr,
		int need_shuffle,
		int *errflag)
{
    MPI_Status status;
    int        rank, comm_size;
    int mask, src, dst;
    int        mpi_errno = MPI_SUCCESS;
	char *send_start, *recv_start;
    MPI_Comm comm;
	int logp;
	MPI_Status sta;
    MPIU_CHKLMEM_DECL(1);
	

	int nbytes;
	int j = 0;
    
    comm = comm_ptr->handle;
    comm_size = comm_ptr->local_size;
    rank = comm_ptr->rank;
	
	logp = get_logp(comm_size);
	mask = comm_size/2;

	dst = get_bwr(rank, logp);
	
	if(need_shuffle == 1 && rank != dst)
	{
		send_start = recvbuf + rank * recvbytes;
		MPIC_Sendrecv(sendbuf, sendbytes, MPI_BYTE, dst, MPIR_GATHER_TAG, send_start, recvbytes, MPI_BYTE, dst, MPIR_GATHER_TAG, comm_ptr, &sta, errflag);
	}
	else
	{
		send_start = sendbuf;
	}

	nbytes = sendbytes;
	j = logp;

	/* Use binomial-tree algorithm */
	while(mask > 0)
	{
		if((rank>>j) == 0)
		{
			dst = rank ^ mask;
			if(rank < dst)
			{
				recv_start = (char *)(send_start + nbytes);
				MPIC_Recv(recv_start, nbytes, MPI_BYTE, dst, MPIR_GATHER_TAG, comm_ptr, &sta, errflag);
			}
			else
			{
				MPIC_Send(send_start, nbytes, MPI_BYTE, dst, MPIR_GATHER_TAG, comm_ptr, errflag);
			}
		}
		mask >>= 1;
		nbytes <<= 1;
		j--;
	}

fn_exit:
    return mpi_errno;
 fn_fail:
    goto fn_exit;

}


#undef FUNCNAME
#define FUNCNAME MPIC_Allgather_positive
#undef FCNAME
#define FCNAME MPIU_QUOTE(FUNCNAME)
int MPIC_Allgather_positive ( 
		const void *sendbuf,
		int sendbytes,
		void *recvbuf,
		int recvbytes,
		MPID_Comm *comm_ptr,
		int need_shuffle,
		int shuffle_bytes,
		int *errflag )
{
    MPI_Status status;
    int        rank, comm_size;
    int mask, src, dst;
    int        mpi_errno = MPI_SUCCESS;
	char *send_start, *recv_start;
    MPI_Comm comm;
	int logp;
	MPI_Status sta;
    MPIU_CHKLMEM_DECL(1);

	int nbytes;
	int j = 0;
    
    comm = comm_ptr->handle;
    comm_size = comm_ptr->local_size;
    rank = comm_ptr->rank;
	
	logp = get_logp(comm_size);
	mask = 1;

	dst = get_bwr(rank, logp);
	
	send_start = (char *)(recvbuf + rank * recvbytes);

	if(need_shuffle == 1 && rank != dst)
	{
		MPIC_Sendrecv(sendbuf, sendbytes, MPI_BYTE, dst, MPIR_ALLGATHER_TAG, send_start, recvbytes, MPI_BYTE, dst, MPIR_ALLGATHER_TAG, comm_ptr, &sta, errflag);
	}
	else
	{
		memcpy(send_start, sendbuf, sendbytes);
	}

	nbytes = sendbytes;

	/* Use recursive-doubling algorithm */
	while(mask < comm_size)
	{
		dst = rank ^ mask;
		send_start = (char *)(recvbuf + ((rank>>j)<<j) * recvbytes);
		recv_start = (char *)(recvbuf + ((dst>>j)<<j) * sendbytes);//my recv_start is remote send start
		MPIC_Sendrecv(send_start, nbytes, MPI_BYTE, dst, MPIR_ALLGATHER_TAG, recv_start, nbytes, MPI_BYTE, dst, MPIR_ALLGATHER_TAG, comm_ptr, &sta, errflag);
		
		mask <<= 1;
		nbytes <<= 1;
		j++;
	}

fn_exit:
    return mpi_errno;
 fn_fail:
    goto fn_exit;
	
	
}


#undef FUNCNAME
#define FUNCNAME MPIC_Allgather_negative
#undef FCNAME
#define FCNAME MPIU_QUOTE(FUNCNAME)
int MPIC_Allgather_negative ( 
		const void *sendbuf,
		int sendbytes,
		void *recvbuf,
		int recvbytes,
		MPID_Comm *comm_ptr,
		int need_shuffle,
		int shuffle_bytes,
		int *errflag )
{
    MPI_Status status;
    int        rank, comm_size;
    int mask, src, dst;
    int        mpi_errno = MPI_SUCCESS;
	char *send_start, *recv_start;
    MPI_Comm comm;
	int logp;
	MPI_Status sta;
    MPIU_CHKLMEM_DECL(1);

	int nbytes;
    
    comm = comm_ptr->handle;
    comm_size = comm_ptr->local_size;
    rank = comm_ptr->rank;
	
	logp = get_logp(comm_size);
	mask = comm_size/2;

	dst = get_bwr(rank, logp);
	
	send_start = (char *)(recvbuf + dst * recvbytes);

	if(need_shuffle == 1 && rank != dst)
	{
		MPIC_Sendrecv(sendbuf, sendbytes, MPI_BYTE, dst, MPIR_ALLGATHER_TAG, send_start, recvbytes, MPI_BYTE, dst, MPIR_ALLGATHER_TAG, comm_ptr, &sta, errflag);
	}
	else
	{
		memcpy(send_start, sendbuf, sendbytes);
	}

	nbytes = sendbytes;

	/* Use recursive-doubling algorithm */
	while(mask > 0)
	{
		dst = rank ^ mask;
		send_start = (char *)(recvbuf + get_bwr(rank%(mask<<1), logp) * recvbytes);
		recv_start = (char *)(recvbuf + get_bwr(dst%(mask<<1), logp) * sendbytes);//my recv_start is remote send start
		
		double t1, t2;
		MPIC_Sendrecv(send_start, nbytes, MPI_BYTE, dst, MPIR_ALLGATHER_TAG, recv_start, nbytes, MPI_BYTE, dst, MPIR_ALLGATHER_TAG, comm_ptr, &sta, errflag);
		mask >>= 1;
		nbytes <<= 1;
	}

fn_exit:
    return mpi_errno;
 fn_fail:
    goto fn_exit;
	
}

#define LOOP 1000
#define LOOP_S 10

void wd_micro_bench(size_t len)
{
	athread_enter64();
    MPID_Comm *comm_ptr;
    MPID_Comm_get_ptr( MPI_COMM_WORLD, comm_ptr );
	int errflag;
	double t1, t2;
	double tt;
	double tt_min, tt_max, tt_ave;
	int i;

	int *A;
	int *B;
	int logp;
	void *sendbuf;
	logp = get_logp(comm_ptr->local_size);
	if(comm_ptr->rank == 0)
	{
		printf("-------------------- comm_size = %d ----------------------\n", comm_ptr->local_size);
	}
	
	A = valloc(len);
	B = valloc(comm_ptr->local_size * len);
	sendbuf = B + get_bwr(comm_ptr->rank, logp) * (len/sizeof(int));

	if(A == NULL || B == NULL)
	{
		E_PRINT("malloc failed!");
	}
	tt = 0;
	for(i = 0; i < LOOP + LOOP_S; i++)
	{
		t1 = MPI_Wtime();
		MPIC_Gather_positive(sendbuf, len, B, len, 0, comm_ptr, 0, &errflag);
		t2 = MPI_Wtime();
		if(i >= LOOP_S)
		{
			tt += (t2 - t1);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	MPI_Reduce(&tt, &tt_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&tt, &tt_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&tt, &tt_ave, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if(comm_ptr->rank == 0)
	{
		tt_min *= (1000*1000);
		tt_max *= (1000*1000);
		tt_ave *= (1000*1000);
		tt_min /= LOOP;
		tt_max /= LOOP;
		tt_ave /= (LOOP * comm_ptr->local_size);
		printf("Gather_positive:\t\t len = \t%8d KB, max = \t%8.2f us, min = \t%8.2f us, ave = \t%8.2f us\n", (len)/1024, tt_max, tt_min, tt_ave);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	tt = 0;
	for(i = 0; i < LOOP + LOOP_S; i++)
	{
		t1 = MPI_Wtime();
		MPIC_Gather_positive(sendbuf, len, B, len, 0, comm_ptr, 1, &errflag);
		t2 = MPI_Wtime();
		if(i >= LOOP_S)
		{
			tt += (t2 - t1);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	MPI_Reduce(&tt, &tt_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&tt, &tt_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&tt, &tt_ave, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if(comm_ptr->rank == 0)
	{
		tt_min *= (1000*1000);
		tt_max *= (1000*1000);
		tt_ave *= (1000*1000);
		tt_min /= LOOP;
		tt_max /= LOOP;
		tt_ave /= (LOOP * comm_ptr->local_size);
		printf("Gather_positive(shuffle):\t len = \t%8d KB, max = \t%8.2f us, min = \t%8.2f us, ave = \t%8.2f us\n", (len)/1024, tt_max, tt_min, tt_ave);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	tt = 0;
	for(i = 0; i < LOOP + LOOP_S; i++)
	{
		t1 = MPI_Wtime();
		MPIC_Gather_negative(sendbuf, len, B, len, 0, comm_ptr, 0, &errflag);
		t2 = MPI_Wtime();
		if(i >= LOOP_S)
		{
			tt += (t2 - t1);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	MPI_Reduce(&tt, &tt_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&tt, &tt_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&tt, &tt_ave, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if(comm_ptr->rank == 0)
	{
		tt_min *= (1000*1000);
		tt_max *= (1000*1000);
		tt_ave *= (1000*1000);
		tt_min /= LOOP;
		tt_max /= LOOP;
		tt_ave /= (LOOP * comm_ptr->local_size);
		printf("Gather_negative:\t\t len = \t%8d KB, max = \t%8.2f us, min = \t%8.2f us, ave = \t%8.2f us\n", (len)/1024, tt_max, tt_min, tt_ave);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	tt = 0;
	for(i = 0; i < LOOP + LOOP_S; i++)
	{
		t1 = MPI_Wtime();
		MPIC_Allgather_positive(sendbuf, len, B, len, comm_ptr, 0, 0, &errflag);
		t2 = MPI_Wtime();
		if(i >= LOOP_S)
		{
			tt += (t2 - t1);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	MPI_Reduce(&tt, &tt_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&tt, &tt_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&tt, &tt_ave, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if(comm_ptr->rank == 0)
	{
		tt_min *= (1000*1000);
		tt_max *= (1000*1000);
		tt_ave *= (1000*1000);
		tt_min /= LOOP;
		tt_max /= LOOP;
		tt_ave /= (LOOP * comm_ptr->local_size);
		printf("Allgather_positive:\t\t len = \t%8d KB, max = \t%8.2f us, min = \t%8.2f us, ave = \t%8.2f us\n", (len)/1024, tt_max, tt_min, tt_ave);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	tt = 0;
	for(i = 0; i < LOOP + LOOP_S; i++)
	{
		t1 = MPI_Wtime();
		MPIC_Allgather_negative(sendbuf, len, B, len, comm_ptr, 0, 0, &errflag);
		t2 = MPI_Wtime();
		if(i >= LOOP_S)
		{
			tt += (t2 - t1);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	MPI_Reduce(&tt, &tt_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&tt, &tt_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&tt, &tt_ave, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if(comm_ptr->rank == 0)
	{
		tt_min *= (1000*1000);
		tt_max *= (1000*1000);
		tt_ave *= (1000*1000);
		tt_min /= LOOP;
		tt_max /= LOOP;
		tt_ave /= (LOOP * comm_ptr->local_size);
		printf("Allgather_negative:\t\t len = \t%8d KB, max = \t%8.2f us, min = \t%8.2f us, ave = \t%8.2f us\n", (len)/1024, tt_max, tt_min, tt_ave);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	tt = 0;
	for(i = 0; i < LOOP + LOOP_S; i++)
	{
		t1 = MPI_Wtime();
		MPIC_Allgather_negative(sendbuf, len, B, len, comm_ptr, 1, 0, &errflag);
		t2 = MPI_Wtime();
		if(i >= LOOP_S)
		{
			tt += (t2 - t1);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	MPI_Reduce(&tt, &tt_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&tt, &tt_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&tt, &tt_ave, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if(comm_ptr->rank == 0)
	{
		tt_min *= (1000*1000);
		tt_max *= (1000*1000);
		tt_ave *= (1000*1000);
		tt_min /= LOOP;
		tt_max /= LOOP;
		tt_ave /= (LOOP * comm_ptr->local_size);
		printf("Allgather_negative(Shuffle):\t len = \t%8d KB, max = \t%8.2f us, min = \t%8.2f us, ave = \t%8.2f us\n", (len)/1024, tt_max, tt_min, tt_ave);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(comm_ptr->rank == 0)
	{
		printf("\n\n\n");fflush(stdout);
	}
	
	free(A);
	free(B);

	athread_leave64();
}
