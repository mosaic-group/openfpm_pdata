#!groovy

timeout(180)
{

parallel (


"nyu" : {node ('nyu')
                  {
                    deleteDir()

                    int ntry = 5
                    while (ntry != 0)
                    {
                      try {
                        checkout scm
                        ntry = 0
                      }
                      catch (IOException e)
                      {
                        ntry--
                        sleep(50)
                      }
                    }


                    stage ('build_nyu')
                    {
                      sh "./build.sh $WORKSPACE $NODE_NAME pdata 0 $BRANCH_NAME }"
                    }

                    stage ('run_nyu')
                    {
                      parallel (
                      "1" : {sh "./run.sh $WORKSPACE $NODE_NAME 1 0 0  $BRANCH_NAME"},
                      "2" : {sh "./run.sh $WORKSPACE $NODE_NAME 2 0 0  $BRANCH_NAME"},
                      "3" : {sh "./run.sh $WORKSPACE $NODE_NAME 3 0 0  $BRANCH_NAME"})
                      sh "./run.sh $WORKSPACE $NODE_NAME 5 0 0 $BRANCH_NAME"
                      sh "./success.sh 2 nyu opefpm_pdata"
                    }
                  }
                 },




"sb15" : {node ('sbalzarini-mac-15')
                  {
                    deleteDir()
                    env.PATH = "/usr/local/bin:${env.PATH}"
               
                    int ntry = 5
                    while (ntry != 0)
                    {
                      try {
                        checkout scm
                        ntry = 0
                      }
                      catch (IOException e)
                      {
                        ntry--
                        sleep(50)
                      }
                    }


                    stage ('build_sb15')
                    {
                      sh "./build.sh $WORKSPACE $NODE_NAME pdata 0  $BRANCH_NAME"
                    }

                    stage ('run_sb15')
                    {
                      parallel (
                      "1" : {sh "./run.sh $WORKSPACE $NODE_NAME 1 0 0  $BRANCH_NAME"},
                      "2" : {sh "./run.sh $WORKSPACE $NODE_NAME 2 0 0  $BRANCH_NAME"},
                      "3" : {sh "./run.sh $WORKSPACE $NODE_NAME 3 0 0  $BRANCH_NAME"}
                      )
                      sh "./run.sh $WORKSPACE $NODE_NAME 4 0 0 $BRANCH_NAME"
                      sh "./run.sh $WORKSPACE $NODE_NAME 5 0 0 $BRANCH_NAME"
                      sh "./run.sh $WORKSPACE $NODE_NAME 6 0 0 $BRANCH_NAME"
                      sh "./run.sh $WORKSPACE $NODE_NAME 7 0 0 $BRANCH_NAME"
                      sh "./success.sh 2 sbalzarini-mac-15 opefpm_pdata"
                    }
                  }
                 },


"gin" : {node ('gin')
                  {
                    deleteDir()
                    env.PATH = "/usr/local/bin:${env.PATH}"

                    int ntry = 5
                    while (ntry != 0)
                    {
                      try {
                        checkout scm
                        ntry = 0
                      }
                      catch (IOException e)
                      {
                        ntry--
                        sleep(50)
                      }
                    }


                    stage ('build_gin')
                    {
                      sh "./build.sh $WORKSPACE $NODE_NAME pdata 0 $BRANCH_NAME"
                    }

                    stage ('run_gin')
                    {
                      parallel (
                      "p1" : {sh "./run.sh $WORKSPACE $NODE_NAME 1 0 0  $BRANCH_NAME"},
                      "p2" : {sh "./run.sh $WORKSPACE $NODE_NAME 2 0 0  $BRANCH_NAME"},
                      "p3" : {sh "./run.sh $WORKSPACE $NODE_NAME 3 0 0  $BRANCH_NAME"},
                      "p4" : {sh "./run.sh $WORKSPACE $NODE_NAME 5 0 0  $BRANCH_NAME"}
                      )
                      sh "./success.sh 2 gin opefpm_pdata"
                    }
                  }
                 }

)
}

