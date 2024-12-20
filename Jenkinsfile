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
					dockerImage = docker.build("${DOCKER_IMG}", "--pull .")
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
		stage('DeployProdStack') {
			agent {
				label 'helm'
			}
	        when {
                buildingTag();
	        }
		    steps {
	            script {
                    withKubeConfig([credentialsId: "prod-kube-config"]) {
                        sh('helm repo update && helm dependency update helm-chart/banzai-nres/ '+
                                '&& helm package helm-chart/banzai-nres --app-version="${GIT_DESCRIPTION}" --version="${GIT_DESCRIPTION}" ' +
                                '&& helm upgrade --install banzai-nres banzai-nres-"${GIT_DESCRIPTION}".tgz --namespace=prod ' +
                                '--set image.tag="${GIT_DESCRIPTION}" --values=helm-chart/banzai-nres/values-prod.yaml ' +
                                '--force --atomic --timeout=3600s')
                    }
                 }
		    }
		}
		stage('DeployTestStack') {
			agent {
				label 'helm'
			}
			when {
				anyOf {
					branch 'PR-*'
					branch 'dev'
					expression { return params.forceEndToEnd }
				}
			}
		    steps {
	            script {
                    withKubeConfig([credentialsId: "build-kube-config"]) {
                        if (env.BRANCH_NAME == "dev") {
                            dataTag = '1.0.4'
                        } else {
                            dataTag = '1.0.4-slim'
                        }
                        sh('helm repo update')
                        final cmd = " helm delete --namespace=build --purge banzai-nres-e2e &> cleanup.txt"
                        final status = sh(script: cmd, returnStatus: true)
                        final output = readFile('cleanup.txt').trim()
                        sh(script: "rm -f cleanup.txt", returnStatus: true)
                        echo output
                        sh(script: "kubectl delete pvc banzai-nres-e2e --wait=true --timeout=600s", returnStatus: true)
                        sh('helm upgrade --namespace=build --install banzai-nres-e2e helm-chart/banzai-nres-e2e ' +
                            '--set banzaiNRES.tag="${GIT_DESCRIPTION}" --set dataImage.tag=' + dataTag +
                            ' --force --wait --timeout=3600s')

                        podName = sh(script: 'kubectl get po -l app.kubernetes.io/instance=banzai-nres-e2e ' +
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
                    withKubeConfig([credentialsId: "build-kube-config"]) {
						sh("kubectl exec ${podName} -c banzai-nres-e2e-listener -- " +
						        "pytest -s --durations=0 --junitxml=/home/archive/pytest-master-bias.xml " +
						        "-m master_bias /lco/banzai-nres/")
					}
				}
			}
			post {
				always {
					script {
					    withKubeConfig([credentialsId: "build-kube-config"]) {
						    sh("kubectl cp -c banzai-nres-e2e-listener ${podName}:/home/archive/pytest-master-bias.xml " +
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
                    withKubeConfig([credentialsId: "build-kube-config"]) {
						sh("kubectl exec ${podName} -c banzai-nres-e2e-listener -- " +
						        "pytest -s --durations=0 --junitxml=/home/archive/pytest-master-dark.xml " +
						        "-m master_dark /lco/banzai-nres/")
					}
				}
			}
			post {
				always {
					script {
					    withKubeConfig([credentialsId: "build-kube-config"]) {
						    sh("kubectl cp -c banzai-nres-e2e-listener ${podName}:/home/archive/pytest-master-dark.xml " +
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
                    withKubeConfig([credentialsId: "build-kube-config"]) {
						sh("kubectl exec ${podName} -c banzai-nres-e2e-listener -- " +
						        "pytest -s --durations=0 --junitxml=/home/archive/pytest-master-flat.xml " +
						        "-m master_flat /lco/banzai-nres/")
					}
				}
			}
			post {
				always {
					script {
					    withKubeConfig([credentialsId: "build-kube-config"]) {
						    sh("kubectl cp -c banzai-nres-e2e-listener ${podName}:/home/archive/pytest-master-flat.xml " +
						            "pytest-master-flat.xml")
						    junit "pytest-master-flat.xml"
						}
					}
				}
			}
		}
		stage('Test-Master-Arc-Creation') {
			when {
				anyOf {
					branch 'PR-*'
					expression { return params.forceEndToEnd }
				}
			}
			steps {
				script {
                    withKubeConfig([credentialsId: "build-kube-config"]) {
						sh("kubectl exec ${podName} -c banzai-nres-e2e-listener -- " +
						        "pytest -s --durations=0 --junitxml=/home/archive/pytest-master-arc.xml " +
						        "-m master_arc /lco/banzai-nres/")
					}
				}
			}
			post {
				always {
					script {
					    withKubeConfig([credentialsId: "build-kube-config"]) {
						    sh("kubectl cp -c banzai-nres-e2e-listener ${podName}:/home/archive/pytest-master-arc.xml " +
						            "pytest-master-arc.xml")
						    junit "pytest-master-arc.xml"
						}
					}
				}
			}
		}
		stage('Test-Science-Frame-Creation') {
						agent {
				label 'helm'
			}
			when {
				anyOf {
					branch 'PR-*'
					expression { return params.forceEndToEnd }
				}
			}
			steps {
				script {
                    withKubeConfig([credentialsId: "build-kube-config"]) {
						sh("kubectl exec ${podName} -c banzai-nres-e2e-listener -- " +
						        "pytest -s --durations=0 --junitxml=/home/archive/pytest-science-frames.xml " +
						        "-m science_frames /lco/banzai-nres/")
					}
				}
			}
			post {
				always {
					script {
					    withKubeConfig([credentialsId: "build-kube-config"]) {
						    sh("kubectl cp -c banzai-nres-e2e-listener ${podName}:/home/archive/pytest-science-frames.xml " +
						            "pytest-science-frames.xml")
						    junit "pytest-science-frames.xml"
						}
					}
				}
				success {
					script {
					    withKubeConfig([credentialsId: "build-kube-config"]) {
                            sh("helm delete --namespace=build banzai-nres-e2e --purge || true")
					    }
					}
				}
			}
		}
	}
}
