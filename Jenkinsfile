#!groovy

parallel (


"nyu" : {node ('nyu')
                  {
                    deleteDir()
                    checkout scm
                    stage ('build_nyu')
                    {
                      sh "./build.sh $WORKSPACE $NODE_NAME pdata env.BRANCH_NAME"
                    }

                    stage ('run_nyu')
                    {
                      parallel (
                      "1" : {sh "./run.sh $WORKSPACE $NODE_NAME 1"},
                      "2" : {sh "./run.sh $WORKSPACE $NODE_NAME 2"},
                      "3" : {sh "./run.sh $WORKSPACE $NODE_NAME 3"})
                      sh "./run.sh $WORKSPACE $NODE_NAME 5"
                      sh "./success.sh 2 nyu opefpm_pdata"
                    }
                  }
                 },




"sb15" : {node ('sbalzarini-mac-15')
                  {
                    deleteDir()
                    env.PATH = "/usr/local/bin:${env.PATH}"
                    checkout scm
                    stage ('build_sb15')
                    {
                      sh "./build.sh $WORKSPACE $NODE_NAME pdata env.BRANCH_NAME"
                    }

                    stage ('run_sb15')
                    {
                      parallel (
                      "1" : {sh "./run.sh $WORKSPACE $NODE_NAME 1"},
                      "2" : {sh "./run.sh $WORKSPACE $NODE_NAME 2"},
                      "3" : {sh "./run.sh $WORKSPACE $NODE_NAME 3"}
                      )
                      sh "./run.sh $WORKSPACE $NODE_NAME 4"
                      sh "./run.sh $WORKSPACE $NODE_NAME 5"
                      sh "./run.sh $WORKSPACE $NODE_NAME 6"
                      sh "./run.sh $WORKSPACE $NODE_NAME 7"
                      sh "./success.sh 2 sbalzarini-mac-15 opefpm_pdata"
                    }
                  }
                 },


"gin" : {node ('gin')
                  {
                    deleteDir()
                    env.PATH = "/usr/local/bin:${env.PATH}"
                    checkout scm
                    stage ('build_gin')
                    {
                      sh "./build.sh $WORKSPACE $NODE_NAME pdata env.BRANCH_NAME"
                    }

                    stage ('run_gin')
                    {
                      parallel (
                      "p1" : {sh "./run.sh $WORKSPACE $NODE_NAME 1"},
                      "p2" : {sh "./run.sh $WORKSPACE $NODE_NAME 2"},
                      "p3" : {sh "./run.sh $WORKSPACE $NODE_NAME 3"},
                      "p4" : {sh "./run.sh $WORKSPACE $NODE_NAME 5"}
                      )
                      sh "./success.sh 2 gin opefpm_pdata"
                    }
                  }
                 }

)

