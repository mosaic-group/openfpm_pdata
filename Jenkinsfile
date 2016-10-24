#!groovy

parallel (


"nyu" : {node ('nyu')
                  {
                    deleteDir()
                    checkout scm
                    stage ('build_nyu')
                    {
                      sh "./build.sh $WORKSPACE $NODE_NAME pdata"
                    }

                    stage ('run_nyu')
                    {
                      sh "./run.sh $WORKSPACE $NODE_NAME 1"
                      sh "./run.sh $WORKSPACE $NODE_NAME 2"
                      sh "./run.sh $WORKSPACE $NODE_NAME 3"
                      sh "./run.sh $WORKSPACE $NODE_NAME 5"
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
                      sh "./build.sh $WORKSPACE $NODE_NAME pdata"
                    }

                    stage ('run_sb15')
                    {
                      sh "./run.sh $WORKSPACE $NODE_NAME 1"
                      sh "./run.sh $WORKSPACE $NODE_NAME 2"
                      sh "./run.sh $WORKSPACE $NODE_NAME 3"
                      sh "./run.sh $WORKSPACE $NODE_NAME 4"
                      sh "./run.sh $WORKSPACE $NODE_NAME 5"
                      sh "./run.sh $WORKSPACE $NODE_NAME 6"
                      sh "./run.sh $WORKSPACE $NODE_NAME 7"
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
                      sh "./build.sh $WORKSPACE $NODE_NAME pdata"
                    }

                    stage ('run_gin')
                    {
                      parallel (
                      "p1" : {sh "./run.sh $WORKSPACE $NODE_NAME 1"},
                      "p2" : {sh "./run.sh $WORKSPACE $NODE_NAME 2"},
                      "p3" : {sh "./run.sh $WORKSPACE $NODE_NAME 3"},
                      "p4" : {sh "./run.sh $WORKSPACE $NODE_NAME 5"}
                      )
                    }
                  }
                 }

)

