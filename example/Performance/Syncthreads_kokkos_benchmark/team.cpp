#include <Kokkos_Core.hpp>
#include <cstdio>



int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);

  printf("LayoutRight\n");

  {



    //////////////////////// TEAMS ////////////////////////

    using Kokkos::TeamPolicy;
    using Kokkos::parallel_for;

    typedef TeamPolicy<Kokkos::OpenMP>::member_type member_type;
    // Create an instance of the policy
    int team_sz = 1;
    int sz = 512;
    TeamPolicy<Kokkos::OpenMP> policy (sz*sz, team_sz);
    // Launch a kernel
    
    Kokkos::fence();
    Kokkos::Timer timer;  

    parallel_for (policy, KOKKOS_LAMBDA (member_type team_member) {
        // Calculate a global thread id
         int k = team_member.league_rank () * team_member.team_size () +
                team_member.team_rank ();
        // Calculate the sum of the global thread ids of this team
	team_member.team_barrier();
	team_member.team_barrier();
	team_member.team_barrier();
	team_member.team_barrier();
	team_member.team_barrier();
	team_member.team_barrier();

        team_member.team_barrier();
        team_member.team_barrier();
        team_member.team_barrier();
        team_member.team_barrier();
        team_member.team_barrier();
        team_member.team_barrier();

        team_member.team_barrier();
        team_member.team_barrier();
        team_member.team_barrier();
        team_member.team_barrier();
        team_member.team_barrier();
        team_member.team_barrier();

        team_member.team_barrier();
        team_member.team_barrier();
        team_member.team_barrier();
        team_member.team_barrier();
        team_member.team_barrier();
        team_member.team_barrier();

         // Atomically add the value to a global value
      });

      Kokkos::fence();
        double time = timer.seconds();
        std::cout << "TIME: " << time / (sz*sz*team_sz*24) * 1e9  << " ns"  << std::endl; 

    ///////////////////////////////////////////////////////

  }
  printf("LayoutLeft\n");

  Kokkos::finalize();
}

