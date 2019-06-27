#!/usr/bin/env groovy

@Library('lco-shared-libs@0.0.8') _

pipeline {
	agent any
	parameters {
		booleanParam(
			name: 'forceEndToEnd',
			defaultValue: false,
			description: 'When true, forces the end-to-end tests to always run.')
	}
	environment {
		dockerImage = null
		PROJ_NAME = projName()
		GIT_DESCRIPTION = gitDescription()
		DOCKER_IMG = dockerImageName("${LCO_DOCK_REG}", "${PROJ_NAME}", "${GIT_DESCRIPTION}")
	}
	options {
		timeout(time: 3, unit: 'HOURS')
	    lock resource: 'BANZAINRESLock'
	}
	stages {
		stage('Build image') {
			steps {
				script {
					dockerImage = docker.build("${DOCKER_IMG}")
				}
			}
		}
		stage('Push image') {
			steps {
				script {
					dockerImage.push("${GIT_DESCRIPTION}")
				}
			}
		}
		stage('Unit Tests') {
			steps {
				script {
					sh 'docker run --rm --user=root --entrypoint=pytest ${DOCKER_IMG} -m "not e2e" /lco/banzai-nres/'
				}
			}
		}
		stage('DeployTestStack') {
			when {
				anyOf {
					branch 'PR-*'
					expression { return params.forceEndToEnd }
				}
			}
		    steps {
	            script {
                    withKubeConfig([credentialsId: "dev-kube-config"]) {
                        sh('helm repo update && helm upgrade --install banzai-nres lco/banzai-nres ' +
                                '--set banzaiNRES.tag="${GIT_DESCRIPTION}" --namespace dev --wait --timeout=3600')

                        podName = sh(script: 'kubectl -n dev get po -l app.kubernetes.io/instance=banzai-nres ' +
                                        '--sort-by=.status.startTime -o jsonpath="{.items[-1].metadata.name}"',
                                     returnStdout: true).trim()

                    }
                 }
		    }
		}
		stage('Test-Master-Bias-Creation') {
			when {
				anyOf {
					branch 'PR-*'
					expression { return params.forceEndToEnd }
				}
			}
			steps {
				script {
                    withKubeConfig([credentialsId: "dev-kube-config"]) {
						sh("kubectl exec ${podName} -c banzai-nres-listener -- " +
						        "pytest --durations=0 --junitxml=/home/archive/pytest-master-bias.xml " +
						        "-m master_bias /lco/banzai-nres/")
					}
				}
			}
			post {
				always {
					script {
					    withKubeConfig([credentialsId: "dev-kube-config"]) {
						    sh("kubectl cp ${podName} -c " +
						            "banzai-nres-listener:/home/archive/pytest-master-bias.xml " +
						            "pytest-master-bias.xml")
						    junit "pytest-master-bias.xml"
						}
					}
				}
			}
		}
		stage('Test-Master-Dark-Creation') {
			when {
				anyOf {
					branch 'PR-*'
					expression { return params.forceEndToEnd }
				}
			}
			steps {
				script {
                    withKubeConfig([credentialsId: "dev-kube-config"]) {
						sh("kubectl exec ${podName} -c banzai-nres-listener  -- " +
						        "pytest --durations=0 " --junitxml=/home/archive/pytest-master-dark.xml " +
						        "-m master_dark /lco/banzai-nres/")
					}
				}
			}
			post {
				always {
					script {
					    withKubeConfig([credentialsId: "dev-kube-config"]) {
						    sh("kubectl cp ${podName} -c " +
						            "banzai-nres-listener:/home/archive/pytest-master-dark.xml " +
						            "pytest-master-dark.xml")
						    junit "pytest-master-dark.xml"
						}
					}
				}
			}
		}
		stage('Test-Master-Flat-Creation') {
			when {
				anyOf {
					branch 'PR-*'
					expression { return params.forceEndToEnd }
				}
			}
			steps {
				script {
                    withKubeConfig([credentialsId: "dev-kube-config"]) {
						sh("kubectl exec ${podName} -c banzai-nres -- " +
						        "pytest --durations=0 --junitxml=/home/archive/pytest-master-flat.xml " +
						        "-m master_flat /lco/banzai-nres/")
					}
				}
			}
			post {
				always {
					script {
					    withKubeConfig([credentialsId: "dev-kube-config"]) {
						    sh("kubectl cp ${podName} -c " +
						            "banzai-nres-listener:/home/archive/pytest-master-flat.xml " +
						            "pytest-master-flat.xml")
						    junit "pytest-master-flat.xml"
						}
					}
				}
				success {
					script {
					    withKubeConfig([credentialsId: "dev-kube-config"]) {
                            sh("helm delete banzai-nres || true")
					    }
					}
				}
			}
		}
	}
}
